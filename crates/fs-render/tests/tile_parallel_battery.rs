//! Deterministic tile-render integration battery (bead
//! `frankensim-h7xu5.5.1`). Each failure names the run, worker/tile policy,
//! pixel, channel, and exact binary64 bits needed for replay.

#![cfg(feature = "tracer")]

use asupersync::types::Budget;
use core::mem::size_of;
use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::NumericalCertificate;
use fs_exec::{CancelGate, Cx, ExecMode, RunError, RunId, StreamKey};
use fs_geom::fixtures::SphereChart;
use fs_geom::{Aabb, Chart, ChartSample, Point3, TraceStepClaim, Vec3};
use fs_render::lighting::EnvironmentMap;
use fs_render::spectral::lift_rgb;
use fs_render::tracer::{
    Camera, DirectStrategy, Film, FilmTimeMode, Material, PendingRender, Primitive,
    RenderExecutionConfig, RenderExecutionError, RenderShardError, RenderShardLimits,
    RenderShardMergeLimits, RenderTileLayout, RenderWorkerPool, Sampler, Scene, Settings, Shape,
    TracerError, UniformRenderShardResult, UniformRenderShardSpec, merge_uniform_shards, render,
    render_range_with_execution, render_static_shard, render_with_execution,
};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

const SEED: u64 = 0x7469_6c65_5f76_3101;
const MEMORY_LIMIT_BYTES: u64 = 64 << 20;

fn with_cx<R>(cancelled: bool, operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_cx_mode(cancelled, ExecMode::Deterministic, operation)
}

fn with_cx_mode<R>(cancelled: bool, mode: ExecMode, operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_cx_budget_mode(cancelled, mode, Budget::INFINITE, operation)
}

fn with_cx_budget_mode<R>(
    cancelled: bool,
    mode: ExecMode,
    budget: Budget,
    operation: impl FnOnce(&Cx<'_>) -> R,
) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: SEED,
                kernel_id: 0x5449_4c45,
                tile: 0,
                iteration: 0,
            },
            budget,
            mode,
        );
        operation(&cx)
    })
}

fn scene() -> Scene {
    Scene {
        primitives: Vec::new(),
        lights: Vec::new(),
        environment: Some(
            EnvironmentMap::try_from_linear_srgb(
                4,
                2,
                vec![
                    [0.3, 0.4, 0.5],
                    [0.7, 0.2, 0.1],
                    [0.2, 0.8, 0.4],
                    [0.9, 0.7, 0.3],
                    [0.4, 0.1, 0.8],
                    [0.6, 0.6, 0.6],
                    [0.1, 0.3, 0.9],
                    [0.8, 0.5, 0.2],
                ],
                0.37,
            )
            .expect("valid deterministic environment"),
        ),
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 0.0),
            forward: Vec3::new(1.0, 0.0, 0.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.7,
        },
    }
}

fn settings(spp: u32, sampler: Sampler) -> Settings {
    Settings {
        width: 35,
        height: 19,
        spp,
        max_depth: 2,
        sampler,
        strategy: DirectStrategy::Mis,
        seed: SEED,
    }
}

fn execution(tile_width: u32, tile_height: u32, workers: usize, run: u64) -> RenderExecutionConfig {
    RenderExecutionConfig::try_new(
        tile_width,
        tile_height,
        workers,
        MEMORY_LIMIT_BYTES,
        RunId(run),
    )
    .expect("valid execution policy")
}

fn assert_film_bits_eq(left: &Film, right: &Film, context: &str) {
    assert_eq!(
        (left.width, left.height),
        (right.width, right.height),
        "{context}"
    );
    assert_eq!(left.spp_done, right.spp_done, "{context}");
    assert_eq!(left.time_mode, right.time_mode, "{context}");
    for (pixel, (left, right)) in left.xyz.iter().zip(&right.xyz).enumerate() {
        for channel in 0..3 {
            assert_eq!(
                left[channel].to_bits(),
                right[channel].to_bits(),
                "{context}: pixel={pixel} channel={channel} left={:#018x} right={:#018x}",
                left[channel].to_bits(),
                right[channel].to_bits()
            );
        }
    }
}

struct CancellingSphere {
    evaluations: Arc<AtomicUsize>,
    cancel_at: Option<usize>,
    gate: Option<Arc<CancelGate>>,
}

impl Chart for CancellingSphere {
    fn eval(&self, point: Point3, cx: &Cx<'_>) -> ChartSample {
        let evaluation = self.evaluations.fetch_add(1, Ordering::SeqCst) + 1;
        if self.cancel_at == Some(evaluation)
            && let Some(gate) = &self.gate
        {
            gate.request();
        }
        SphereChart {
            center: Point3::new(3.0, 0.0, 0.0),
            radius: 1.0,
        }
        .eval(point, cx)
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(2.0, -1.0, -1.0), Point3::new(4.0, 1.0, 1.0))
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::ExactDistance
    }

    fn name(&self) -> &'static str {
        "tile-battery-cancelling-sphere"
    }
}

struct PanickingChart;

impl Chart for PanickingChart {
    fn eval(&self, _point: Point3, _cx: &Cx<'_>) -> ChartSample {
        panic!("declared renderer tile panic")
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(2.0, -1.0, -1.0), Point3::new(4.0, 1.0, 1.0))
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
        "tile-battery-panicking-chart"
    }
}

struct StreamCheckingSphere {
    expected_seed: u64,
    saw_wrong_seed: Arc<AtomicBool>,
}

impl Chart for StreamCheckingSphere {
    fn eval(&self, point: Point3, cx: &Cx<'_>) -> ChartSample {
        if cx.stream_key().seed != self.expected_seed {
            self.saw_wrong_seed.store(true, Ordering::SeqCst);
        }
        SphereChart {
            center: Point3::new(3.0, 0.0, 0.0),
            radius: 1.0,
        }
        .eval(point, cx)
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(2.0, -1.0, -1.0), Point3::new(4.0, 1.0, 1.0))
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::ExactDistance
    }

    fn name(&self) -> &'static str {
        "tile-battery-stream-checking-sphere"
    }
}

struct PanicOnceSphere {
    evaluations: Arc<AtomicUsize>,
    panic_at: usize,
}

impl Chart for PanicOnceSphere {
    fn eval(&self, point: Point3, cx: &Cx<'_>) -> ChartSample {
        let evaluation = self.evaluations.fetch_add(1, Ordering::SeqCst) + 1;
        assert_ne!(
            evaluation, self.panic_at,
            "declared one-time pending-render panic"
        );
        SphereChart {
            center: Point3::new(3.0, 0.0, 0.0),
            radius: 1.0,
        }
        .eval(point, cx)
    }

    fn support(&self) -> Aabb {
        Aabb::new(Point3::new(2.0, -1.0, -1.0), Point3::new(4.0, 1.0, 1.0))
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::ExactDistance
    }

    fn name(&self) -> &'static str {
        "tile-battery-panic-once-sphere"
    }
}

fn with_chart(mut scene: Scene, chart: Box<dyn Chart>) -> Scene {
    scene.primitives.push(Primitive {
        shape: Shape::Chart(chart),
        material: Material::Lambertian {
            reflectance: lift_rgb([0.6, 0.6, 0.6]),
        },
        emission: None,
    });
    scene
}

#[test]
fn g0_tile_layout_covers_odd_image_exactly() {
    let layout = RenderTileLayout::try_new(35, 19, 16, 8).expect("valid odd layout");
    assert_eq!(
        (layout.tiles_x(), layout.tiles_y(), layout.tile_count()),
        (3, 3, 9)
    );
    assert_eq!(
        layout.bounds(8).expect("last tile"),
        fs_render::tracer::RenderTileBounds {
            x: 32,
            y: 16,
            width: 3,
            height: 3,
        }
    );
    let mut visits = vec![0u8; 35 * 19];
    for tile in 0..layout.tile_count() {
        let bounds = layout.bounds(tile).expect("planned tile");
        for y in bounds.y..bounds.y + bounds.height {
            for x in bounds.x..bounds.x + bounds.width {
                visits[y as usize * 35 + x as usize] += 1;
            }
        }
    }
    assert!(visits.iter().all(|visits| *visits == 1));
    assert!(layout.bounds(9).is_none());
    assert!(RenderTileLayout::try_new(u32::MAX, 2, 16, 8).is_err());
    assert!(RenderTileLayout::try_new(1, 1, 0, 1).is_err());
    assert_eq!(
        Film::try_new(u32::MAX, 2),
        Err(TracerError::InvalidInput),
        "Film must reject dimensions outside tracer-v1's u32 pixel identity domain"
    );
}

#[test]
fn g0_one_tile_reports_only_the_worker_actually_admitted() {
    let scene = scene();
    let settings = settings(1, Sampler::Iid);
    let policy = execution(4_096, 4_096, 8, 0x0a11);
    let output = with_cx(false, |cx| {
        render_with_execution(&scene, cx, &settings, &policy)
    })
    .expect("single-tile render");
    assert_eq!(output.report.requested_workers, 8);
    assert_eq!(output.report.workers, 1);
    assert_eq!(output.report.executor.tiles_by_worker.len(), 1);
    assert!(output.report.idle_worker_ns <= output.report.traversal_ns);
}

#[test]
fn g0_tile_scratch_envelope_is_admitted_before_dispatch() {
    let scene = scene();
    let settings = settings(1, Sampler::Iid);
    let film_bytes =
        u64::from(settings.width) * u64::from(settings.height) * size_of::<[f64; 3]>() as u64;
    let scratch_bytes = 4 * 8 * 8 * size_of::<[f64; 3]>() as u64;
    let policy =
        RenderExecutionConfig::try_new(8, 8, 4, film_bytes + scratch_bytes - 1, RunId(0x0a12))
            .expect("valid near-limit policy");
    let result = with_cx(false, |cx| {
        render_with_execution(&scene, cx, &settings, &policy)
    });
    let Err(RenderExecutionError::Memory(refusal)) = result else {
        panic!("scratch envelope should refuse before executor dispatch: {result:?}");
    };
    assert_eq!(refusal.what, "render-tile-scratch-envelope");
    assert_eq!(refusal.requested_bytes, scratch_bytes);
    assert_eq!(refusal.used_bytes, film_bytes);
}

#[test]
fn g5_parked_crew_reuses_workers_across_render_jobs() {
    let scene = scene();
    let settings = settings(3, Sampler::Iid);
    let serial = with_cx(false, |cx| render(&scene, cx, &settings)).expect("serial oracle");
    let initial = execution(8, 5, 4, 0x7000);
    let pool = RenderWorkerPool::new(&initial, ExecMode::Deterministic, SEED);
    pool.with_parked_crew_local(|renderer| {
        for ordinal in 0..3 {
            let job = execution(8, 5, 4, 0x7000 + ordinal);
            let output = with_cx(false, |cx| renderer.render(&scene, cx, &settings, &job))
                .expect("parked render job");
            assert_film_bits_eq(
                &serial,
                &output.film,
                &format!("parked render ordinal={ordinal}"),
            );
            assert_eq!(output.report.executor.declared_run, job.run_id());
            assert_eq!(output.report.executor.completed, 20);
            assert_eq!(output.report.tile_scratch_envelope_bytes, 4 * 8 * 5 * 24);
        }

        let explicit_equal_weights = execution(8, 5, 4, 0x70f0)
            .with_quantum_weights(vec![1; 4])
            .expect("explicit equal weights");
        let output = with_cx(false, |cx| {
            renderer.render(&scene, cx, &settings, &explicit_equal_weights)
        })
        .expect("empty and explicit all-one worker weights are equivalent");
        assert_film_bits_eq(&serial, &output.film, "canonical equal worker weights");

        let mismatched = execution(8, 5, 2, 0x7fff);
        let error = with_cx(false, |cx| {
            renderer.render(&scene, cx, &settings, &mismatched)
        })
        .expect_err("worker mismatch must refuse before dispatch");
        assert!(matches!(
            error,
            RenderExecutionError::Config(
                fs_render::tracer::RenderExecutionConfigError::ParkedCrewMismatch
            )
        ));
    });
}

#[test]
fn g5_parked_scheduler_seed_is_not_visible_to_scene_charts() {
    let saw_wrong_seed = Arc::new(AtomicBool::new(false));
    let scene = with_chart(
        scene(),
        Box::new(StreamCheckingSphere {
            expected_seed: SEED,
            saw_wrong_seed: Arc::clone(&saw_wrong_seed),
        }),
    );
    let settings = settings(2, Sampler::Iid);
    let policy = execution(8, 5, 4, 0x7050);
    let serial = with_cx(false, |cx| render(&scene, cx, &settings)).expect("serial oracle");
    let pool = RenderWorkerPool::new(&policy, ExecMode::Deterministic, 0xdead_beef_cafe_f00d);
    let tiled = pool
        .with_parked_crew_local(|renderer| {
            with_cx(false, |cx| renderer.render(&scene, cx, &settings, &policy))
        })
        .expect("parked render with a deliberately unrelated scheduler seed");
    assert_film_bits_eq(&serial, &tiled.film, "parked scheduler seed changed film");
    assert!(
        !saw_wrong_seed.load(Ordering::SeqCst),
        "a scene chart observed the pool placement seed instead of Settings::seed"
    );
}

#[test]
fn g5_pending_render_binds_mode_budget_and_counts_refused_attempts() {
    let scene = scene();
    let settings = settings(1, Sampler::Iid);
    let policy = execution(8, 5, 2, 0x7051);
    let pending = with_cx(false, |cx| {
        PendingRender::begin_static(&scene, cx, settings, policy.clone())
    })
    .expect("deterministic pending render");

    let first = with_cx_mode(false, ExecMode::Fast, |cx| pending.resume(cx))
        .expect_err("mode-changing retry must refuse before dispatch");
    assert!(matches!(
        first.cause(),
        RenderExecutionError::Config(
            fs_render::tracer::RenderExecutionConfigError::ResumeModeMismatch {
                expected: ExecMode::Deterministic,
                actual: ExecMode::Fast,
            }
        )
    ));
    assert_eq!(first.progress().attempts, 1);
    assert_eq!(first.attempt_report().attempt_index, 1);
    assert_eq!(first.attempt_report().executor.total, 0);

    let second = with_cx_mode(false, ExecMode::Fast, |cx| first.into_pending().resume(cx))
        .expect_err("a second mode-changing retry must also refuse");
    assert_eq!(second.progress().attempts, 2);
    assert_eq!(second.attempt_report().attempt_index, 2);
    assert_eq!(second.attempt_report().executor.total, 0);

    let third = with_cx_budget_mode(
        false,
        ExecMode::Deterministic,
        Budget::new().with_cost_quota(65_536),
        |cx| second.into_pending().resume(cx),
    )
    .expect_err("budget-changing retry must refuse before dispatch");
    assert!(matches!(
        third.cause(),
        RenderExecutionError::Config(
            fs_render::tracer::RenderExecutionConfigError::ResumeBudgetMismatch
        )
    ));
    assert_eq!(third.progress().attempts, 3);
    assert_eq!(third.attempt_report().attempt_index, 3);
    assert_eq!(third.attempt_report().executor.total, 0);

    let output = with_cx(false, |cx| third.into_pending().resume(cx))
        .expect("the originally bound mode must remain resumable");
    assert_eq!(output.report.attempt_index, 4);
    assert_eq!(output.report.executor.declared_run, policy.run_id());
}

#[test]
fn g5_parallel_workers_and_schedules_are_bit_exact_to_serial() {
    let scene = scene();
    for sampler in [Sampler::Iid, Sampler::OwenSobol] {
        let settings = settings(5, sampler);
        let serial = with_cx(false, |cx| render(&scene, cx, &settings)).expect("serial oracle");
        for workers in [1, 2, 4, 8] {
            let uniform = execution(16, 8, workers, workers as u64);
            let output = with_cx(false, |cx| {
                render_with_execution(&scene, cx, &settings, &uniform)
            })
            .expect("parallel render");
            assert_film_bits_eq(
                &serial,
                &output.film,
                &format!("sampler={sampler:?} workers={workers} tile=16x8"),
            );
            assert_eq!(output.report.executor.completed, 9);
            assert_eq!(output.report.executor.total, 9);
            assert_eq!(output.report.memory.used_bytes, 0);
            let film_bytes = u64::from(settings.width)
                * u64::from(settings.height)
                * size_of::<[f64; 3]>() as u64;
            assert_eq!(output.report.retained_film_bytes, 0);
            assert_eq!(output.report.staging_film_bytes, film_bytes);
            assert_eq!(
                output.report.sampler_state_bytes,
                if sampler == Sampler::OwenSobol {
                    3 * size_of::<[u32; 32]>() as u64
                } else {
                    0
                },
                "only Owen-Sobol jobs may allocate direction state"
            );
            if workers == 1 {
                assert!(
                    output.report.memory.peak_bytes < 2 * film_bytes,
                    "single-worker fresh render retained a second film: peak={} film_bytes={film_bytes}",
                    output.report.memory.peak_bytes
                );
            }
            println!(
                "{{\"suite\":\"fs-render/tile-parallel\",\"case\":\"serial-equality\",\"sampler\":\"{sampler:?}\",\"workers\":{workers},\"tiles\":{},\"setup_ns\":{},\"traversal_ns\":{},\"compute_ns\":{},\"merge_ns\":{},\"idle_worker_ns\":{},\"memory_peak_bytes\":{}}}",
                output.report.layout.tile_count(),
                output.report.setup_ns,
                output.report.traversal_ns,
                output.report.tile_compute_ns,
                output.report.tile_merge_ns,
                output.report.idle_worker_ns,
                output.report.memory.peak_bytes,
            );

            let weights = (0..workers).map(|worker| (worker + 1) as u32).collect();
            let skewed = execution(7, 5, workers, 0x100 + workers as u64)
                .with_quantum_weights(weights)
                .expect("valid skewed weights");
            let output = with_cx(false, |cx| {
                render_with_execution(&scene, cx, &settings, &skewed)
            })
            .expect("skewed parallel render");
            assert_film_bits_eq(
                &serial,
                &output.film,
                &format!("sampler={sampler:?} workers={workers} tile=7x5 skewed"),
            );
        }
    }
}

#[test]
fn g5_progressive_partitions_can_change_worker_and_tile_policy() {
    let scene = scene();
    let settings = settings(7, Sampler::Iid);
    let serial = with_cx(false, |cx| render(&scene, cx, &settings)).expect("serial oracle");
    let mut film = Film::new(settings.width, settings.height);
    for (range_index, (from, to, tile_width, tile_height, workers)) in
        [(0, 2, 16, 8, 4), (2, 5, 7, 5, 2), (5, 7, 32, 32, 8)]
            .into_iter()
            .enumerate()
    {
        let policy = execution(tile_width, tile_height, workers, u64::from(from));
        let report = with_cx(false, |cx| {
            render_range_with_execution(&scene, cx, &settings, &mut film, from, to, &policy)
        })
        .expect("progressive tile range");
        assert_eq!(report.memory.used_bytes, 0);
        if range_index == 0 {
            let film_bytes = u64::from(settings.width)
                * u64::from(settings.height)
                * size_of::<[f64; 3]>() as u64;
            assert!(
                report.memory.peak_bytes >= 2 * film_bytes,
                "progressive transaction failed to account retained+staging films"
            );
            assert_eq!(report.retained_film_bytes, film_bytes);
            assert_eq!(report.staging_film_bytes, film_bytes);
        }
    }
    assert_film_bits_eq(&serial, &film, "parallel progressive 2+3+2 != serial 7");
}

#[test]
fn g4_precancel_and_memory_refusal_leave_film_unchanged() {
    let scene = scene();
    let settings = settings(3, Sampler::Iid);
    let mut film = Film::new(settings.width, settings.height);
    let before = film.clone();
    let policy = execution(8, 8, 4, 44);
    let cancelled = with_cx(true, |cx| {
        render_range_with_execution(&scene, cx, &settings, &mut film, 0, 3, &policy)
    });
    assert!(matches!(
        cancelled,
        Err(fs_render::tracer::RenderExecutionError::Tracer(
            fs_render::tracer::TracerError::Cancelled
        ))
    ));
    assert_film_bits_eq(&before, &film, "pre-cancelled render published pixels");

    let refused = RenderExecutionConfig::try_new(8, 8, 4, 1, RunId(45))
        .expect("one-byte policy is structurally valid");
    let result = with_cx(false, |cx| {
        render_range_with_execution(&scene, cx, &settings, &mut film, 0, 3, &refused)
    });
    assert!(matches!(
        result,
        Err(fs_render::tracer::RenderExecutionError::Memory(_))
    ));
    assert_film_bits_eq(&before, &film, "memory refusal published pixels");
}

#[test]
fn g4_mid_trace_cancel_drains_then_retry_replays_exactly() {
    let settings = settings(3, Sampler::Iid);
    let policy = execution(8, 5, 4, 50);
    let gate = Arc::new(CancelGate::new_clock_free());
    let evaluations = Arc::new(AtomicUsize::new(0));
    let cancelling_scene = with_chart(
        scene(),
        Box::new(CancellingSphere {
            evaluations: Arc::clone(&evaluations),
            cancel_at: Some(1),
            gate: Some(Arc::clone(&gate)),
        }),
    );
    let mut film = Film::new(settings.width, settings.height);
    let before = film.clone();
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let result = arenas.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: SEED,
                kernel_id: 0x4341_4e43,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        render_range_with_execution(&cancelling_scene, &cx, &settings, &mut film, 0, 3, &policy)
    });
    assert!(matches!(
        result,
        Err(RenderExecutionError::Tracer(TracerError::Cancelled))
    ));
    assert!(evaluations.load(Ordering::SeqCst) >= 1);
    assert_film_bits_eq(&before, &film, "mid-trace cancellation published a tile");

    with_cx(false, |cx| {
        render_range_with_execution(&cancelling_scene, cx, &settings, &mut film, 0, 3, &policy)
    })
    .expect("retry under a fresh cancellation authority");
    let reference_scene = with_chart(
        scene(),
        Box::new(CancellingSphere {
            evaluations: Arc::new(AtomicUsize::new(0)),
            cancel_at: None,
            gate: None,
        }),
    );
    let reference = with_cx(false, |cx| render(&reference_scene, cx, &settings))
        .expect("serial retry reference");
    assert_film_bits_eq(&reference, &film, "cancelled retry drifted from serial");
}

#[test]
fn g4_owned_row_prefix_resume_uses_one_film_and_never_double_counts() {
    let settings = settings(3, Sampler::Iid);
    let policy = execution(8, 5, 4, 0x5155);
    let gate = Arc::new(CancelGate::new_clock_free());
    let evaluations = Arc::new(AtomicUsize::new(0));
    let resumable_scene = with_chart(
        scene(),
        Box::new(CancellingSphere {
            evaluations: Arc::clone(&evaluations),
            cancel_at: Some(2_000),
            gate: Some(Arc::clone(&gate)),
        }),
    );
    let pending = with_cx(false, |cx| {
        PendingRender::begin_static(&resumable_scene, cx, settings, policy.clone())
    })
    .expect("admit owned pending render");
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let suspended = arenas
        .scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: SEED,
                    kernel_id: 0x5253_554d,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            pending.resume(&cx)
        })
        .expect_err("injected cancellation must suspend the private job");
    assert!(matches!(
        suspended.cause(),
        RenderExecutionError::Tracer(TracerError::Cancelled)
    ));
    assert!(suspended.attempt_report().layout.tile_count() > 1);
    assert_eq!(suspended.attempt_report().requested_workers, 4);
    assert_eq!(suspended.attempt_report().workers, 4);
    let partial = suspended.progress();
    assert!(
        partial.committed_tile_rows > 0 && partial.committed_tile_rows < partial.total_tile_rows,
        "cancellation must retain a strict row prefix for this fixture: {partial:?} evaluations={}",
        evaluations.load(Ordering::SeqCst)
    );
    assert_eq!(
        suspended.attempt_report().memory.used_bytes,
        suspended.attempt_report().staging_film_bytes
            + suspended.attempt_report().progress_state_bytes
            + suspended.attempt_report().sampler_state_bytes
    );
    assert_eq!(
        suspended.attempt_report().executor.declared_run,
        policy.run_id()
    );

    let retry_pool = RenderWorkerPool::new(&policy, ExecMode::Deterministic, SEED);
    let output = retry_pool
        .with_parked_crew_local(|renderer| {
            with_cx(false, |cx| {
                suspended.into_pending().resume_on_parked(renderer, cx)
            })
        })
        .expect("fresh-authority parked resume must finish");
    let reference_scene = with_chart(
        scene(),
        Box::new(CancellingSphere {
            evaluations: Arc::new(AtomicUsize::new(0)),
            cancel_at: None,
            gate: None,
        }),
    );
    let reference = with_cx(false, |cx| render(&reference_scene, cx, &settings))
        .expect("serial owned-resume reference");
    assert_film_bits_eq(&reference, &output.film, "owned row-prefix resume drifted");
    assert_eq!(output.report.attempt_index, 2);
    assert_eq!(output.report.executor.declared_run, policy.run_id());
    assert_eq!(output.report.memory.used_bytes, 0);
    assert_eq!(output.report.retained_film_bytes, 0);
    let film_bytes =
        u64::from(settings.width) * u64::from(settings.height) * size_of::<[f64; 3]>() as u64;
    assert_eq!(output.report.staging_film_bytes, film_bytes);
    assert!(
        output.report.memory.peak_bytes < 2 * film_bytes,
        "owned resume exceeded one-film plus bounded metadata/scratch: peak={} film={film_bytes}",
        output.report.memory.peak_bytes
    );
}

#[test]
fn g4_owned_row_prefix_survives_contained_worker_panic() {
    let settings = settings(3, Sampler::Iid);
    let policy = execution(35, 19, 1, 0x5158);
    let evaluations = Arc::new(AtomicUsize::new(0));
    let resumable_scene = with_chart(
        scene(),
        Box::new(PanicOnceSphere {
            evaluations: Arc::clone(&evaluations),
            panic_at: 2_000,
        }),
    );
    let pending = with_cx(false, |cx| {
        PendingRender::begin_static(&resumable_scene, cx, settings, policy.clone())
    })
    .expect("admit panic-resumable render");
    let suspended = with_cx(false, |cx| pending.resume(cx))
        .expect_err("the declared one-time worker panic must suspend the job");
    assert!(matches!(
        suspended.cause(),
        RenderExecutionError::Executor(RunError::TilePanicked { .. })
    ));
    assert_eq!(
        suspended.attempt_report().executor.declared_run,
        policy.run_id()
    );
    let partial = suspended.progress();
    assert!(
        partial.committed_tile_rows > 0 && partial.committed_tile_rows < partial.total_tile_rows,
        "panic fixture must leave a strict committed prefix: {partial:?} evaluations={}",
        evaluations.load(Ordering::SeqCst)
    );

    let output = with_cx(false, |cx| suspended.into_pending().resume(cx))
        .expect("fresh attempt must resume after the contained one-time panic");
    let reference_scene = with_chart(
        scene(),
        Box::new(CancellingSphere {
            evaluations: Arc::new(AtomicUsize::new(0)),
            cancel_at: None,
            gate: None,
        }),
    );
    let reference = with_cx(false, |cx| render(&reference_scene, cx, &settings))
        .expect("serial panic-resume reference");
    assert_film_bits_eq(
        &reference,
        &output.film,
        "panic resume double-counted a row",
    );
    assert_eq!(output.report.attempt_index, 2);
    assert_eq!(output.report.executor.declared_run, policy.run_id());
    assert_eq!(output.report.memory.used_bytes, 0);
}

#[test]
fn g0_owned_resume_scratch_refusal_dispatches_no_tiles_and_retains_zero_progress() {
    let scene = scene();
    let settings = settings(1, Sampler::Iid);
    let layout = RenderTileLayout::try_new(settings.width, settings.height, 8, 8)
        .expect("owned refusal layout");
    let film_bytes =
        u64::from(settings.width) * u64::from(settings.height) * size_of::<[f64; 3]>() as u64;
    let persistent_bytes = film_bytes + layout.tile_count() * size_of::<u32>() as u64;
    let policy = RenderExecutionConfig::try_new(8, 8, 4, persistent_bytes, RunId(0x5156))
        .expect("persistent-only policy");
    let pending = with_cx(false, |cx| {
        PendingRender::begin_static(&scene, cx, settings, policy)
    })
    .expect("persistent job state fits exactly");
    let suspended = with_cx(false, |cx| pending.resume(cx))
        .expect_err("row scratch must refuse before executor dispatch");
    let progress = suspended.progress();
    assert_eq!(progress.committed_tile_rows, 0);
    assert_eq!(progress.completed_tiles, 0);
    assert_eq!(suspended.attempt_report().executor.total, 0);
    assert!(matches!(suspended.cause(), RenderExecutionError::Memory(_)));
}

#[test]
fn g4_precompleted_zero_sample_job_survives_cancel_before_publication() {
    let scene = scene();
    let settings = settings(0, Sampler::OwenSobol);
    let policy = execution(8, 8, 2, 0x5157);
    let pending = with_cx(false, |cx| {
        PendingRender::begin_static(&scene, cx, settings, policy)
    })
    .expect("zero-sample pending job");
    assert_eq!(pending.progress().completed_tiles, 15);
    let suspended = with_cx(true, |cx| pending.resume(cx))
        .expect_err("final publication must honor the cancelled authority");
    assert_eq!(suspended.progress().completed_tiles, 15);
    assert_eq!(suspended.attempt_report().executor.total, 0);
    let output = with_cx(false, |cx| suspended.into_pending().resume(cx))
        .expect("fresh authority publishes without retracing");
    assert_eq!(output.report.attempt_index, 2);
    assert_eq!(output.report.executor.total, 0);
    assert_eq!(output.report.sampler_state_bytes, 0);
    assert_eq!(output.film.spp_done, 0);
    assert!(output.film.xyz.iter().all(|pixel| *pixel == [0.0; 3]));
}

#[test]
fn g4_worker_panic_is_contained_and_film_remains_transactional() {
    let settings = settings(2, Sampler::Iid);
    let policy = execution(8, 5, 4, 60);
    let panicking_scene = with_chart(scene(), Box::new(PanickingChart));
    let mut film = Film::new(settings.width, settings.height);
    let before = film.clone();
    let result = with_cx(false, |cx| {
        render_range_with_execution(&panicking_scene, cx, &settings, &mut film, 0, 2, &policy)
    });
    assert!(matches!(
        result,
        Err(RenderExecutionError::Executor(
            RunError::TilePanicked { .. }
        ))
    ));
    assert_film_bits_eq(&before, &film, "contained worker panic published a tile");

    let healthy_scene = scene();
    with_cx(false, |cx| {
        render_range_with_execution(&healthy_scene, cx, &settings, &mut film, 0, 2, &policy)
    })
    .expect("renderer remains reusable after a contained panic");
    let reference = with_cx(false, |cx| render(&healthy_scene, cx, &settings))
        .expect("serial post-panic reference");
    assert_film_bits_eq(&reference, &film, "post-panic retry drifted from serial");
}

fn shard_test_identity(label: &str) -> ContentHash {
    hash_domain("org.frankensim.test.render-sharding.v1", label.as_bytes())
}

fn shard_test_settings(sampler: Sampler) -> Settings {
    Settings {
        width: 6,
        height: 5,
        spp: 4,
        max_depth: 2,
        sampler,
        strategy: DirectStrategy::Mis,
        seed: SEED,
    }
}

#[allow(clippy::too_many_arguments)]
fn shard_spec(
    plan_identity: ContentHash,
    frame_identity: ContentHash,
    settings: Settings,
    layout: RenderTileLayout,
    tile_start: u64,
    tile_end: u64,
    sample_start: u32,
    sample_end: u32,
    limits: RenderShardLimits,
) -> UniformRenderShardSpec {
    UniformRenderShardSpec::try_new(
        plan_identity,
        frame_identity,
        17,
        settings,
        FilmTimeMode::Static,
        layout,
        tile_start,
        tile_end,
        sample_start,
        sample_end,
        limits,
    )
    .expect("valid bounded static render shard")
}

fn render_shard_set(
    scene: &Scene,
    specs: &[UniformRenderShardSpec],
    cx: &Cx<'_>,
) -> Vec<UniformRenderShardResult> {
    specs
        .iter()
        .map(|spec| render_static_shard(scene, cx, spec).expect("deterministic shard render"))
        .collect()
}

#[test]
fn g0_uniform_shards_enforce_caps_codec_completeness_and_transactional_retry() {
    let scene = scene();
    let settings = shard_test_settings(Sampler::Iid);
    let layout = RenderTileLayout::try_new(settings.width, settings.height, 3, 3)
        .expect("four-tile irregular-edge layout");
    assert_eq!(layout.tile_count(), 4);
    let plan_identity = shard_test_identity("bounded-plan");
    let frame_identity = shard_test_identity("bounded-frame");
    let generous = RenderShardLimits::try_new(1 << 20, 4 << 20).expect("generous shard caps");

    let probe = shard_spec(
        plan_identity,
        frame_identity,
        settings,
        layout,
        0,
        3,
        0,
        2,
        generous,
    );
    assert_eq!(probe.payload_pixel_count(), 24);
    assert_eq!(probe.path_count(), 48);
    let exact_limits = RenderShardLimits::try_new(probe.path_count(), probe.encoded_result_bytes())
        .expect("positive exact shard caps");
    let exact_spec = shard_spec(
        plan_identity,
        frame_identity,
        settings,
        layout,
        0,
        3,
        0,
        2,
        exact_limits,
    );
    assert_eq!(exact_spec.path_count(), exact_limits.max_paths());
    assert_eq!(
        exact_spec.encoded_result_bytes(),
        exact_limits.max_result_bytes()
    );
    assert!(matches!(
        UniformRenderShardSpec::try_new(
            plan_identity,
            frame_identity,
            17,
            settings,
            FilmTimeMode::Static,
            layout,
            0,
            3,
            0,
            2,
            RenderShardLimits::try_new(probe.path_count() - 1, 4 << 20)
                .expect("one-short path policy"),
        ),
        Err(RenderShardError::PathLimit { .. })
    ));
    assert!(matches!(
        UniformRenderShardSpec::try_new(
            plan_identity,
            frame_identity,
            17,
            settings,
            FilmTimeMode::Static,
            layout,
            0,
            3,
            0,
            2,
            RenderShardLimits::try_new(1 << 20, probe.encoded_result_bytes() - 1)
                .expect("one-short result-byte policy"),
        ),
        Err(RenderShardError::ResultByteLimit { .. })
    ));

    let exact_result = with_cx(false, |cx| {
        render_static_shard(&scene, cx, &exact_spec).expect("exact-cap shard render")
    });
    let encoded = with_cx(false, |cx| {
        exact_result
            .encode_canonical(exact_result.encoded_result_bytes(), cx)
            .expect("exact-cap canonical shard encoding")
    });
    assert_eq!(encoded.len() as u64, exact_result.encoded_result_bytes());
    assert!(
        with_cx(false, |cx| exact_result
            .encode_canonical(exact_result.encoded_result_bytes() - 1, cx))
        .is_err()
    );
    let decoded = with_cx(false, |cx| {
        UniformRenderShardResult::decode_canonical(
            &encoded,
            encoded.len() as u64,
            &exact_spec,
            exact_spec.plan_identity(),
            exact_spec.shard_identity(),
            cx,
        )
        .expect("strict pinned shard decode")
    });
    assert_eq!(decoded, exact_result);
    for prefix in [0, 1, 7, encoded.len() - 1] {
        assert!(
            with_cx(false, |cx| UniformRenderShardResult::decode_canonical(
                &encoded[..prefix],
                encoded.len() as u64,
                &exact_spec,
                exact_spec.plan_identity(),
                exact_spec.shard_identity(),
                cx,
            ))
            .is_err(),
            "accepted truncated shard prefix {prefix}/{}",
            encoded.len()
        );
    }
    let mut corrupt = encoded.clone();
    let corrupt_index = corrupt.len() / 2;
    corrupt[corrupt_index] ^= 1;
    assert!(
        with_cx(false, |cx| UniformRenderShardResult::decode_canonical(
            &corrupt,
            corrupt.len() as u64,
            &exact_spec,
            exact_spec.plan_identity(),
            exact_spec.shard_identity(),
            cx,
        ))
        .is_err()
    );
    let mut trailing = encoded.clone();
    trailing.push(0);
    assert!(matches!(
        with_cx(false, |cx| UniformRenderShardResult::decode_canonical(
            &trailing,
            trailing.len() as u64,
            &exact_spec,
            exact_spec.plan_identity(),
            exact_spec.shard_identity(),
            cx,
        )),
        Err(RenderShardError::ResultByteLimit { .. } | RenderShardError::TrailingBytes)
    ));
    assert!(matches!(
        with_cx(false, |cx| UniformRenderShardResult::decode_canonical(
            &encoded,
            encoded.len() as u64,
            &exact_spec,
            shard_test_identity("wrong-plan-pin"),
            exact_spec.shard_identity(),
            cx,
        )),
        Err(RenderShardError::PlanIdentityMismatch { .. })
    ));
    assert!(matches!(
        with_cx(false, |cx| UniformRenderShardResult::decode_canonical(
            &encoded,
            encoded.len() as u64,
            &exact_spec,
            exact_spec.plan_identity(),
            shard_test_identity("wrong-shard-pin"),
            cx,
        )),
        Err(RenderShardError::ShardIdentityMismatch { .. })
    ));

    let specs = [
        shard_spec(
            plan_identity,
            frame_identity,
            settings,
            layout,
            0,
            2,
            0,
            2,
            generous,
        ),
        shard_spec(
            plan_identity,
            frame_identity,
            settings,
            layout,
            0,
            2,
            2,
            4,
            generous,
        ),
        shard_spec(
            plan_identity,
            frame_identity,
            settings,
            layout,
            2,
            4,
            0,
            2,
            generous,
        ),
        shard_spec(
            plan_identity,
            frame_identity,
            settings,
            layout,
            2,
            4,
            2,
            4,
            generous,
        ),
    ];
    let results = with_cx(false, |cx| render_shard_set(&scene, &specs, cx));
    let mut generously_capped_bytes = with_cx(false, |cx| {
        results[0]
            .encode_canonical(generous.max_result_bytes(), cx)
            .expect("generously capped canonical result")
    });
    generously_capped_bytes.push(0);
    assert!(matches!(
        with_cx(false, |cx| UniformRenderShardResult::decode_canonical(
            &generously_capped_bytes,
            generously_capped_bytes.len() as u64,
            &specs[0],
            specs[0].plan_identity(),
            specs[0].shard_identity(),
            cx,
        )),
        Err(RenderShardError::TrailingBytes)
    ));
    let input_bytes = results
        .iter()
        .map(UniformRenderShardResult::encoded_result_bytes)
        .sum::<u64>();
    let output_bytes = u64::from(settings.width) * u64::from(settings.height) * 24;
    let exact_merge_limits =
        RenderShardMergeLimits::try_new(input_bytes, output_bytes).expect("exact merge caps");
    with_cx(false, |cx| {
        merge_uniform_shards(&specs, &results, exact_merge_limits, cx)
    })
    .expect("complete result set at exact aggregate caps");
    let one_short_input = RenderShardMergeLimits::try_new(input_bytes - 1, output_bytes)
        .expect("one-short input cap");
    for submitted in [&results[..], &[][..]] {
        match with_cx(false, |cx| {
            merge_uniform_shards(&specs, submitted, one_short_input, cx)
        }) {
            Err(RenderShardError::AggregateInputLimit { limit, observed }) => {
                assert_eq!(limit, input_bytes - 1);
                assert_eq!(observed, input_bytes);
            }
            other => panic!(
                "one-short exact aggregate cap must precede result inspection, got {other:?}"
            ),
        }
    }
    assert!(matches!(
        with_cx(false, |cx| merge_uniform_shards(
            &specs,
            &results,
            RenderShardMergeLimits::try_new(input_bytes, output_bytes - 1)
                .expect("one-short output cap"),
            cx,
        )),
        Err(RenderShardError::OutputByteLimit { .. })
    ));
    assert!(matches!(
        with_cx(false, |cx| merge_uniform_shards(
            &specs,
            &results[..results.len() - 1],
            RenderShardMergeLimits::try_new(64 << 20, 64 << 20).expect("merge policy"),
            cx,
        )),
        Err(RenderShardError::MissingShard(_))
    ));
    assert!(matches!(
        with_cx(false, |cx| merge_uniform_shards(
            &specs[1..],
            &[],
            RenderShardMergeLimits::try_new(64 << 20, 64 << 20).expect("merge policy"),
            cx,
        )),
        Err(RenderShardError::CoverageGap { .. })
    ));
    let overlapping_specs = [specs[0], specs[0], specs[1], specs[2], specs[3]];
    assert!(matches!(
        with_cx(false, |cx| merge_uniform_shards(
            &overlapping_specs,
            &[],
            RenderShardMergeLimits::try_new(64 << 20, 64 << 20).expect("merge policy"),
            cx,
        )),
        Err(RenderShardError::CoverageOverlap { .. })
    ));

    let mut exact_duplicate_results = results.clone();
    exact_duplicate_results.push(results[0].clone());
    with_cx(false, |cx| {
        merge_uniform_shards(&specs, &exact_duplicate_results, exact_merge_limits, cx)
    })
    .expect("exact duplicate is idempotent and consumes no second aggregate-byte charge");
    let mut altered_scene = crate::scene();
    altered_scene.environment = Some(
        EnvironmentMap::try_from_linear_srgb(4, 2, vec![[0.1, 0.9, 0.2]; 8], 0.11)
            .expect("alternate finite emitter"),
    );
    let conflict = with_cx(false, |cx| {
        render_static_shard(&altered_scene, cx, &specs[0]).expect("alternate valid shard result")
    });
    assert_ne!(conflict.result_identity(), results[0].result_identity());
    let mut conflicting_results = results.clone();
    conflicting_results.push(conflict);
    assert!(matches!(
        with_cx(false, |cx| merge_uniform_shards(
            &specs,
            &conflicting_results,
            RenderShardMergeLimits::try_new(64 << 20, 64 << 20).expect("merge policy"),
            cx,
        )),
        Err(RenderShardError::ConflictingDuplicate(_))
    ));

    let foreign_plan_spec = shard_spec(
        shard_test_identity("foreign-plan"),
        frame_identity,
        settings,
        layout,
        0,
        2,
        0,
        2,
        generous,
    );
    let foreign_plan_result = with_cx(false, |cx| {
        render_static_shard(&scene, cx, &foreign_plan_spec).expect("foreign-plan result")
    });
    let mut foreign_results = results.clone();
    foreign_results.push(foreign_plan_result);
    assert!(matches!(
        with_cx(false, |cx| merge_uniform_shards(
            &specs,
            &foreign_results,
            RenderShardMergeLimits::try_new(64 << 20, 64 << 20).expect("merge policy"),
            cx,
        )),
        Err(RenderShardError::ForeignPlan(_))
    ));
    let foreign_frame_spec = shard_spec(
        plan_identity,
        shard_test_identity("foreign-frame"),
        settings,
        layout,
        0,
        2,
        0,
        2,
        generous,
    );
    let foreign_frame_result = with_cx(false, |cx| {
        render_static_shard(&scene, cx, &foreign_frame_spec).expect("foreign-frame result")
    });
    foreign_results.pop();
    foreign_results.push(foreign_frame_result);
    assert!(matches!(
        with_cx(false, |cx| merge_uniform_shards(
            &specs,
            &foreign_results,
            RenderShardMergeLimits::try_new(64 << 20, 64 << 20).expect("merge policy"),
            cx,
        )),
        Err(RenderShardError::ForeignFrame(_))
    ));
    let unexpected_spec = shard_spec(
        plan_identity,
        frame_identity,
        settings,
        layout,
        0,
        4,
        0,
        settings.spp,
        generous,
    );
    let unexpected_result = with_cx(false, |cx| {
        render_static_shard(&scene, cx, &unexpected_spec).expect("unexpected valid shard result")
    });
    foreign_results.pop();
    foreign_results.push(unexpected_result);
    assert!(matches!(
        with_cx(false, |cx| merge_uniform_shards(
            &specs,
            &foreign_results,
            RenderShardMergeLimits::try_new(64 << 20, 64 << 20).expect("merge policy"),
            cx,
        )),
        Err(RenderShardError::UnexpectedShard(_))
    ));

    assert!(matches!(
        with_cx(true, |cx| render_static_shard(&scene, cx, &specs[0])),
        Err(RenderShardError::Cancelled)
    ));
    let retried = with_cx(false, |cx| {
        render_static_shard(&scene, cx, &specs[0]).expect("fresh-authority retry")
    });
    assert_eq!(retried, results[0]);
}

#[test]
fn g3_uniform_shard_invalid_input_diagnostics_are_permutation_invariant() {
    let scene = scene();
    let settings = shard_test_settings(Sampler::Iid);
    let layout = RenderTileLayout::try_new(settings.width, settings.height, 3, 3)
        .expect("four-tile irregular-edge layout");
    let plan_identity = shard_test_identity("g3-plan");
    let frame_identity = shard_test_identity("g3-frame");
    let limits = RenderShardLimits::try_new(1 << 20, 4 << 20).expect("shard caps");
    let merge_limits = RenderShardMergeLimits::try_new(64 << 20, 64 << 20).expect("merge caps");
    let specs = [(0, 2, 0, 2), (0, 2, 2, 4), (2, 4, 0, 2), (2, 4, 2, 4)].map(
        |(tile_start, tile_end, sample_start, sample_end)| {
            shard_spec(
                plan_identity,
                frame_identity,
                settings,
                layout,
                tile_start,
                tile_end,
                sample_start,
                sample_end,
                limits,
            )
        },
    );

    // The reference spec must come from canonical logical order, not caller
    // order. Reversing this invalid expected set therefore reports the same
    // foreign plan identity.
    let foreign_expected = shard_spec(
        shard_test_identity("g3-foreign-expected-plan"),
        frame_identity,
        settings,
        layout,
        2,
        4,
        2,
        4,
        limits,
    );
    let mut mixed_expected = specs;
    mixed_expected[3] = foreign_expected;
    for expected in [mixed_expected.to_vec(), {
        let mut reversed = mixed_expected.to_vec();
        reversed.reverse();
        reversed
    }] {
        match with_cx(false, |cx| {
            merge_uniform_shards(&expected, &[], merge_limits, cx)
        }) {
            Err(RenderShardError::ForeignPlan(actual)) => {
                assert_eq!(actual, foreign_expected.plan_identity());
            }
            other => panic!("expected-order permutation changed rejection: {other:?}"),
        }
    }

    // Structural result errors use fixed class precedence and the minimum
    // offending authority identity, independent of submission order.
    let foreign_plan_specs =
        ["g3-foreign-result-plan-a", "g3-foreign-result-plan-b"].map(|label| {
            shard_spec(
                shard_test_identity(label),
                frame_identity,
                settings,
                layout,
                0,
                2,
                0,
                2,
                limits,
            )
        });
    let foreign_frame_spec = shard_spec(
        plan_identity,
        shard_test_identity("g3-foreign-result-frame"),
        settings,
        layout,
        0,
        2,
        0,
        2,
        limits,
    );
    let unexpected_spec = shard_spec(
        plan_identity,
        frame_identity,
        settings,
        layout,
        0,
        4,
        0,
        settings.spp,
        limits,
    );
    let mut invalid_results = with_cx(false, |cx| {
        vec![
            render_static_shard(&scene, cx, &foreign_frame_spec).expect("foreign-frame result"),
            render_static_shard(&scene, cx, &foreign_plan_specs[0])
                .expect("first foreign-plan result"),
            render_static_shard(&scene, cx, &unexpected_spec).expect("unexpected result"),
            render_static_shard(&scene, cx, &foreign_plan_specs[1])
                .expect("second foreign-plan result"),
        ]
    });
    let expected_foreign_plan = foreign_plan_specs
        .iter()
        .map(|spec| spec.plan_identity())
        .min()
        .expect("two foreign plan identities");
    for submitted in [invalid_results.clone(), {
        invalid_results.reverse();
        invalid_results
    }] {
        match with_cx(false, |cx| {
            merge_uniform_shards(&specs, &submitted, merge_limits, cx)
        }) {
            Err(RenderShardError::ForeignPlan(actual)) => {
                assert_eq!(actual, expected_foreign_plan);
            }
            other => panic!("result permutation changed structural rejection: {other:?}"),
        }
    }

    // Multiple valid conflicting duplicates likewise report the canonical
    // minimum shard identity rather than whichever conflict arrived first.
    let baseline = with_cx(false, |cx| render_shard_set(&scene, &specs, cx));
    let mut altered_scene = crate::scene();
    altered_scene.environment = Some(
        EnvironmentMap::try_from_linear_srgb(4, 2, vec![[0.8, 0.1, 0.4]; 8], 0.17)
            .expect("alternate finite emitter"),
    );
    let conflicts = with_cx(false, |cx| {
        [
            render_static_shard(&altered_scene, cx, &specs[0]).expect("first alternate result"),
            render_static_shard(&altered_scene, cx, &specs[3]).expect("second alternate result"),
        ]
    });
    assert_ne!(
        conflicts[0].result_identity(),
        baseline[0].result_identity()
    );
    assert_ne!(
        conflicts[1].result_identity(),
        baseline[3].result_identity()
    );
    let expected_conflict = [specs[0].shard_identity(), specs[3].shard_identity()]
        .into_iter()
        .min()
        .expect("two conflicting shard identities");
    let mut conflicting_results = baseline;
    conflicting_results.extend(conflicts);
    for submitted in [conflicting_results.clone(), {
        conflicting_results.reverse();
        conflicting_results
    }] {
        match with_cx(false, |cx| {
            merge_uniform_shards(&specs, &submitted, merge_limits, cx)
        }) {
            Err(RenderShardError::ConflictingDuplicate(actual)) => {
                assert_eq!(actual, expected_conflict);
            }
            other => panic!("result permutation changed conflict rejection: {other:?}"),
        }
    }
}

#[test]
fn g5_uniform_shards_are_arrival_order_invariant_with_explicit_exactness_boundaries() {
    for sampler in [Sampler::Iid, Sampler::OwenSobol] {
        let scene = scene();
        let settings = shard_test_settings(sampler);
        let layout = RenderTileLayout::try_new(settings.width, settings.height, 3, 3)
            .expect("four-tile irregular-edge layout");
        let plan_identity = shard_test_identity(match sampler {
            Sampler::Iid => "iid-plan",
            Sampler::OwenSobol => "sobol-plan",
        });
        let frame_identity = shard_test_identity("g5-frame");
        let limits = RenderShardLimits::try_new(1 << 20, 4 << 20).expect("shard caps");
        let merge_limits = RenderShardMergeLimits::try_new(64 << 20, 64 << 20).expect("merge caps");

        let serial =
            with_cx(false, |cx| render(&scene, cx, &settings)).expect("legacy serial reference");
        let tile_only_specs = [
            shard_spec(
                plan_identity,
                frame_identity,
                settings,
                layout,
                0,
                2,
                0,
                settings.spp,
                limits,
            ),
            shard_spec(
                plan_identity,
                frame_identity,
                settings,
                layout,
                2,
                4,
                0,
                settings.spp,
                limits,
            ),
        ];
        let tile_only_results = with_cx(false, |cx| render_shard_set(&scene, &tile_only_specs, cx));
        let tile_only = with_cx(false, |cx| {
            merge_uniform_shards(&tile_only_specs, &tile_only_results, merge_limits, cx)
        })
        .expect("complete tile-only merge");
        assert_film_bits_eq(
            &serial,
            &tile_only,
            "full-SPP tile-only shards must retain legacy serial bits",
        );

        let split_specs = [
            shard_spec(
                plan_identity,
                frame_identity,
                settings,
                layout,
                0,
                2,
                0,
                2,
                limits,
            ),
            shard_spec(
                plan_identity,
                frame_identity,
                settings,
                layout,
                0,
                2,
                2,
                4,
                limits,
            ),
            shard_spec(
                plan_identity,
                frame_identity,
                settings,
                layout,
                2,
                4,
                0,
                2,
                limits,
            ),
            shard_spec(
                plan_identity,
                frame_identity,
                settings,
                layout,
                2,
                4,
                2,
                4,
                limits,
            ),
        ];
        let split_results = with_cx(false, |cx| render_shard_set(&scene, &split_specs, cx));
        let canonical = with_cx(false, |cx| {
            merge_uniform_shards(&split_specs, &split_results, merge_limits, cx)
        })
        .expect("canonical split-sample merge");
        let mut reversed = split_results.clone();
        reversed.reverse();
        let reverse_merged = with_cx(false, |cx| {
            merge_uniform_shards(&split_specs, &reversed, merge_limits, cx)
        })
        .expect("reverse-arrival split-sample merge");
        assert_film_bits_eq(
            &canonical,
            &reverse_merged,
            "split-sample reverse arrival changed frozen-plan bits",
        );
        let permuted = vec![
            split_results[2].clone(),
            split_results[0].clone(),
            split_results[3].clone(),
            split_results[1].clone(),
            split_results[2].clone(),
        ];
        let permuted_merged = with_cx(false, |cx| {
            merge_uniform_shards(&split_specs, &permuted, merge_limits, cx)
        })
        .expect("permuted split-sample merge with exact duplicate");
        assert_film_bits_eq(
            &canonical,
            &permuted_merged,
            "split-sample permutation or exact duplicate changed frozen-plan bits",
        );
    }

    // Positive regression for an exact cover that is not a Cartesian product
    // of one tile partition and one sample partition. Sample one spans both
    // tiles while sample zero is split at the tile boundary.
    let scene = scene();
    let mut settings = shard_test_settings(Sampler::Iid);
    settings.width = 4;
    settings.height = 1;
    settings.spp = 2;
    let layout = RenderTileLayout::try_new(4, 1, 2, 1).expect("two-tile layout");
    let plan_identity = shard_test_identity("non-cartesian-cover-plan");
    let frame_identity = shard_test_identity("non-cartesian-cover-frame");
    let limits = RenderShardLimits::try_new(1 << 20, 4 << 20).expect("shard caps");
    let specs = [
        shard_spec(
            plan_identity,
            frame_identity,
            settings,
            layout,
            0,
            2,
            1,
            2,
            limits,
        ),
        shard_spec(
            plan_identity,
            frame_identity,
            settings,
            layout,
            1,
            2,
            0,
            1,
            limits,
        ),
        shard_spec(
            plan_identity,
            frame_identity,
            settings,
            layout,
            0,
            1,
            0,
            1,
            limits,
        ),
    ];
    let results = with_cx(false, |cx| render_shard_set(&scene, &specs, cx));
    let merged = with_cx(false, |cx| {
        merge_uniform_shards(
            &specs,
            &results,
            RenderShardMergeLimits::try_new(4 << 20, 4 << 20).expect("merge caps"),
            cx,
        )
    })
    .expect("non-Cartesian exact cover");
    let serial =
        with_cx(false, |cx| render(&scene, cx, &settings)).expect("non-Cartesian serial oracle");
    assert_film_bits_eq(
        &serial,
        &merged,
        "two-tile/two-sample non-Cartesian exact cover changed serial bits",
    );
}
