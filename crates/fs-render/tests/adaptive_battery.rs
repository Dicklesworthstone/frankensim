//! Deterministic adaptive-sampling battery (bead
//! `frankensim-h7xu5.5.2`). Diagnostics include the sampler, execution
//! policy, pixel/channel, and binary64 bits needed to replay any mismatch.

#![cfg(feature = "tracer")]

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, RunId, StreamKey};
use fs_geom::fixtures::SphereChart;
use fs_geom::{Aabb, Chart, ChartSample, Point3, TraceStepClaim, Vec3};
use fs_render::charts::TriMesh;
use fs_render::lighting::EnvironmentMap;
use fs_render::spectral::lift_rgb;
use fs_render::tracer::{
    ADAPTIVE_SAMPLING_SEMANTICS_VERSION, AdaptiveDecision, AdaptiveFilm, AdaptiveSamplingConfig,
    AdaptiveSamplingError, Camera, DirectStrategy, Film, Material, PendingAdaptiveRender,
    Primitive, RectLight, RenderExecutionConfig, RenderExecutionError, RenderTileLayout,
    RenderWorkerPool, Sampler, Scene, Settings, Shape, TracerError, render,
    render_adaptive_with_execution,
};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

const SEED: u64 = 0x6164_6170_745f_7631;
const MEMORY_LIMIT_BYTES: u64 = 64 << 20;

fn with_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    with_gate_cx(&gate, operation)
}

fn with_gate_cx<R>(gate: &CancelGate, operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: SEED,
                kernel_id: 0x4144_4150,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx)
    })
}

fn execution(tile_width: u32, tile_height: u32, workers: usize, run: u64) -> RenderExecutionConfig {
    RenderExecutionConfig::try_new(
        tile_width,
        tile_height,
        workers,
        MEMORY_LIMIT_BYTES,
        RunId(run),
    )
    .expect("valid adaptive execution policy")
}

fn settings(width: u32, height: u32, spp: u32, sampler: Sampler) -> Settings {
    Settings {
        width,
        height,
        spp,
        max_depth: 2,
        sampler,
        strategy: DirectStrategy::Mis,
        seed: SEED,
    }
}

fn camera() -> Camera {
    Camera {
        eye: Point3::new(0.0, 0.0, 0.0),
        forward: Vec3::new(1.0, 0.0, 0.0),
        up: Vec3::new(0.0, 1.0, 0.0),
        half_tan: 0.7,
    }
}

fn black_scene() -> Scene {
    // A physical emitter behind the camera satisfies lighting admission while
    // remaining unreachable by every primary ray in this empty scene. The
    // analytic image is therefore exactly black without a test-only tracer
    // bypass.
    let corner = Point3::new(-10.0, -1.0, -1.0);
    let edge_u = Vec3::new(0.0, 2.0, 0.0);
    let edge_v = Vec3::new(0.0, 0.0, 2.0);
    let emission = (lift_rgb([1.0, 1.0, 1.0]), 1.0);
    let light_mesh = TriMesh::new(
        vec![
            [-10.0, -1.0, -1.0],
            [-10.0, 1.0, -1.0],
            [-10.0, 1.0, 1.0],
            [-10.0, -1.0, 1.0],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    );
    Scene {
        primitives: vec![Primitive {
            shape: Shape::Mesh(light_mesh),
            material: Material::Lambertian {
                reflectance: lift_rgb([0.0, 0.0, 0.0]),
            },
            emission: Some(emission),
        }],
        lights: vec![RectLight {
            corner,
            edge_u,
            edge_v,
            prim: 0,
            emission,
        }],
        environment: None,
        camera: camera(),
    }
}

fn environment_scene() -> Scene {
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
        camera: camera(),
    }
}

fn sparse_glossy_scene() -> Scene {
    let mut scene = black_scene();
    scene.primitives.push(Primitive {
        shape: Shape::Chart(Box::new(SphereChart {
            center: Point3::new(3.0, 0.0, 0.0),
            radius: 0.8,
        })),
        material: Material::Ggx {
            reflectance: lift_rgb([0.9, 0.9, 0.9]),
            alpha: 0.08,
        },
        emission: None,
    });
    scene
}

fn strict_policy(minimum: u32) -> AdaptiveSamplingConfig {
    AdaptiveSamplingConfig::try_new(minimum, 1, 0.0, 0.0, 0.0)
        .expect("valid zero-dispersion-only policy")
}

fn assert_f64x3_bits_eq(left: [f64; 3], right: [f64; 3], context: &str) {
    for channel in 0..3 {
        assert_eq!(
            left[channel].to_bits(),
            right[channel].to_bits(),
            "{context}: channel={channel} left={:#018x} right={:#018x}",
            left[channel].to_bits(),
            right[channel].to_bits()
        );
    }
}

fn assert_adaptive_bits_eq(left: &AdaptiveFilm, right: &AdaptiveFilm, context: &str) {
    assert_eq!(
        (left.width(), left.height()),
        (right.width(), right.height()),
        "{context}"
    );
    assert_eq!(left.maximum_samples(), right.maximum_samples(), "{context}");
    assert_eq!(left.policy(), right.policy(), "{context}");
    assert_eq!(left.sampler(), right.sampler(), "{context}");
    assert_eq!(left.stream_seed(), right.stream_seed(), "{context}");
    assert_eq!(
        left.semantics_version(),
        right.semantics_version(),
        "{context}"
    );
    assert_eq!(left.time_mode(), right.time_mode(), "{context}");
    assert_eq!(left.sample_counts(), right.sample_counts(), "{context}");
    assert_eq!(left.decisions(), right.decisions(), "{context}");
    for pixel in 0..left.xyz_sums().len() {
        assert_f64x3_bits_eq(
            left.xyz_sums()[pixel],
            right.xyz_sums()[pixel],
            &format!("{context}: raw-sum pixel={pixel}"),
        );
        assert_f64x3_bits_eq(
            left.running_means_xyz()[pixel],
            right.running_means_xyz()[pixel],
            &format!("{context}: estimator-mean pixel={pixel}"),
        );
        assert_f64x3_bits_eq(
            left.m2_xyz()[pixel],
            right.m2_xyz()[pixel],
            &format!("{context}: m2 pixel={pixel}"),
        );
    }
}

fn assert_uniform_raw_sum_bits_eq(uniform: &Film, adaptive: &AdaptiveFilm, context: &str) {
    assert_eq!(
        (uniform.width, uniform.height),
        (adaptive.width(), adaptive.height())
    );
    assert_eq!(uniform.xyz.len(), adaptive.xyz_sums().len());
    for (pixel, (uniform, adaptive)) in uniform.xyz.iter().zip(adaptive.xyz_sums()).enumerate() {
        assert_f64x3_bits_eq(*uniform, *adaptive, &format!("{context}: pixel={pixel}"));
    }
}

#[test]
fn g0_constant_scene_stops_at_minimum_with_complete_aovs_and_exact_cost() {
    let scene = black_scene();
    let settings = settings(11, 7, 16, Sampler::Iid);
    let policy = strict_policy(2);
    let execution = execution(4, 3, 3, 0x1001);
    let output =
        with_cx(|cx| render_adaptive_with_execution(&scene, cx, &settings, policy, &execution))
            .expect("constant black adaptive render");
    let film = &output.film;
    let pixels = u64::from(settings.width) * u64::from(settings.height);

    assert_eq!(film.xyz_sums().len(), pixels as usize);
    assert_eq!(film.running_means_xyz().len(), pixels as usize);
    assert_eq!(film.m2_xyz().len(), pixels as usize);
    assert!(film.xyz_sums().iter().all(|xyz| *xyz == [0.0; 3]));
    assert!(film.m2_xyz().iter().all(|xyz| *xyz == [0.0; 3]));
    assert!(film.sample_counts().iter().all(|samples| *samples == 2));
    assert!(
        film.decisions()
            .iter()
            .all(|decision| *decision == AdaptiveDecision::ErrorThreshold)
    );
    assert_eq!(film.maximum_samples(), 16);
    assert_eq!(film.policy(), policy);
    assert_eq!(film.sampler(), Sampler::Iid);
    assert_eq!(film.stream_seed(), SEED);
    assert_eq!(
        film.semantics_version(),
        ADAPTIVE_SAMPLING_SEMANTICS_VERSION
    );
    assert_eq!(film.beauty_mean_xyz(0), Some([0.0; 3]));
    assert_eq!(film.estimator_mean_xyz(0), Some([0.0; 3]));
    assert_eq!(film.sample_variance_xyz(0), Some([0.0; 3]));
    assert_eq!(film.dispersion_proxy_xyz(0), Some([0.0; 3]));
    assert_eq!(film.beauty_mean_xyz(pixels as usize), None);

    let summary = output.summary();
    assert_eq!(summary.pixels, pixels);
    assert_eq!(summary.minimum_samples, 2);
    assert_eq!(summary.maximum_samples, 2);
    assert_eq!(summary.total_samples, pixels * 2);
    assert_eq!(summary.converged_pixels, pixels);
    assert_eq!(summary.maximum_sample_pixels, 0);
    assert_eq!(output.report.memory.used_bytes, 0);

    let layout = RenderTileLayout::try_new(11, 7, 4, 3).unwrap();
    let mut tile_paths = 0_u64;
    for tile in 0..layout.tile_count() {
        let summary = film.tile_summary(layout, tile).expect("owned tile summary");
        assert_eq!(summary.minimum_samples, 2);
        assert_eq!(summary.maximum_samples, 2);
        assert_eq!(summary.total_samples, summary.pixels * 2);
        assert_eq!(summary.maximum_dispersion_xyz, [0.0; 3]);
        tile_paths += summary.total_samples;
    }
    assert_eq!(tile_paths, output.summary().total_samples);
    assert!(
        film.tile_summary(RenderTileLayout::try_new(10, 7, 4, 3).unwrap(), 0)
            .is_none(),
        "a summary layout with different image identity must refuse"
    );

    let uniform = with_cx(|cx| render(&scene, cx, &settings)).expect("uniform black oracle");
    for pixel in 0..pixels as usize {
        assert_f64x3_bits_eq(
            film.beauty_mean_xyz(pixel).unwrap(),
            uniform.xyz[pixel].map(|sum| sum / f64::from(settings.spp)),
            &format!("analytic-black equal-quality oracle pixel={pixel}"),
        );
    }
    assert!(
        output.summary().total_samples < pixels * u64::from(settings.spp),
        "constant analytic reference must use fewer paths than uniform rendering"
    );
}

#[test]
fn g0_invalid_ceiling_and_memory_refusal_happen_before_dispatch() {
    let scene = environment_scene();
    let settings = settings(4, 4, 4, Sampler::Iid);
    let insufficient_memory = RenderExecutionConfig::try_new(2, 2, 2, 1, RunId(0x1002))
        .expect("one-byte limit is structurally valid");
    let refused = with_cx(|cx| {
        render_adaptive_with_execution(
            &scene,
            cx,
            &settings,
            strict_policy(2),
            &insufficient_memory,
        )
    });
    let Err(RenderExecutionError::Memory(refusal)) = refused else {
        panic!("adaptive film allocation should refuse before dispatch: {refused:?}");
    };
    assert_eq!(refusal.what, "render-adaptive-film");
    assert_eq!(refusal.used_bytes, 0);

    let invalid_ceiling = AdaptiveSamplingConfig::try_new(5, 1, 0.0, 0.0, 0.0).unwrap();
    assert_eq!(
        with_cx(|cx| render_adaptive_with_execution(
            &scene,
            cx,
            &settings,
            invalid_ceiling,
            &execution(2, 2, 2, 0x1003),
        )),
        Err(RenderExecutionError::Adaptive(
            AdaptiveSamplingError::MaximumBelowMinimum
        ))
    );
}

#[test]
fn g3_sparse_noisy_fixture_allocates_more_paths_only_where_needed() {
    let scene = sparse_glossy_scene();
    let settings = settings(24, 16, 12, Sampler::Iid);
    let policy = strict_policy(4);
    let output = with_cx(|cx| {
        render_adaptive_with_execution(&scene, cx, &settings, policy, &execution(8, 4, 4, 0x2001))
    })
    .expect("sparse glossy adaptive render");
    let minimum_pixels = output
        .film
        .sample_counts()
        .iter()
        .filter(|samples| **samples == policy.minimum_samples())
        .count();
    let refined_pixels = output
        .film
        .sample_counts()
        .iter()
        .filter(|samples| **samples > policy.minimum_samples())
        .count();
    assert!(
        minimum_pixels > 0,
        "black background never converged at minimum"
    );
    assert!(
        refined_pixels > 0,
        "noisy glossy/silhouette pixels received no extra paths"
    );
    assert!(
        output.summary().total_samples < output.summary().pixels * u64::from(settings.spp),
        "heterogeneous scene did not save any raw paths"
    );

    // Compare against a uniform render whose rounded-up path budget is at
    // least the adaptive cost, then measure both against a disjoint-seed
    // high-SPP reference. This is an empirical fixture-specific result, not a
    // universal threshold or confidence claim.
    let uniform_spp = u32::try_from(
        output
            .summary()
            .total_samples
            .div_ceil(output.summary().pixels),
    )
    .expect("small fixture average spp");
    let uniform_settings = Settings {
        spp: uniform_spp,
        ..settings
    };
    let uniform =
        with_cx(|cx| render(&scene, cx, &uniform_settings)).expect("cost-matched uniform baseline");
    assert!(
        output.summary().total_samples < output.summary().pixels * u64::from(uniform_settings.spp),
        "rounded-up uniform baseline must consume strictly more paths"
    );
    let reference_spp = 32_u32;
    let references = [
        0x9e37_79b9_7f4a_7c15,
        0xd1b5_4a32_d192_ed03,
        0x94d0_49bb_1331_11eb,
    ]
    .map(|seed_mask| {
        let reference_settings = Settings {
            spp: reference_spp,
            seed: SEED ^ seed_mask,
            ..settings
        };
        with_cx(|cx| render(&scene, cx, &reference_settings))
            .expect("independent high-SPP reference replicate")
    });
    let mut adaptive_squared_error = 0.0_f64;
    let mut uniform_squared_error = 0.0_f64;
    for pixel in 0..uniform.xyz.len() {
        let adaptive_mean = output.film.beauty_mean_xyz(pixel).unwrap();
        let uniform_mean = uniform.xyz[pixel].map(|sum| sum / f64::from(uniform_settings.spp));
        let mut reference_mean = [0.0_f64; 3];
        for reference in &references {
            for (mean, sum) in reference_mean.iter_mut().zip(reference.xyz[pixel]) {
                *mean += sum / f64::from(reference_spp);
            }
        }
        reference_mean = reference_mean.map(|sum| sum / 3.0);
        for channel in 0..3 {
            adaptive_squared_error += (adaptive_mean[channel] - reference_mean[channel]).powi(2);
            uniform_squared_error += (uniform_mean[channel] - reference_mean[channel]).powi(2);
        }
    }
    assert!(
        adaptive_squared_error <= uniform_squared_error,
        "adaptive fixture error did not beat the equal-or-higher-cost uniform baseline: adaptive={adaptive_squared_error:?} uniform={uniform_squared_error:?} adaptive_paths={} uniform_paths={}",
        output.summary().total_samples,
        output.summary().pixels * u64::from(uniform_settings.spp),
    );
    println!(
        "{{\"suite\":\"fs-render/adaptive\",\"case\":\"heterogeneous-allocation\",\"sampler\":\"{:?}\",\"seed\":{},\"minimum_spp\":{},\"maximum_spp\":{},\"total_paths\":{},\"uniform_paths\":{},\"minimum_pixels\":{},\"refined_pixels\":{},\"adaptive_reference_squared_error\":{},\"uniform_reference_squared_error\":{}}}",
        settings.sampler,
        settings.seed,
        policy.minimum_samples(),
        settings.spp,
        output.summary().total_samples,
        output.summary().pixels * u64::from(uniform_settings.spp),
        minimum_pixels,
        refined_pixels,
        adaptive_squared_error,
        uniform_squared_error,
    );
}

#[test]
fn g5_adaptive_decisions_and_moments_ignore_workers_tiles_and_parked_reuse() {
    let scene = environment_scene();
    let policy = AdaptiveSamplingConfig::try_new(2, 2, 0.002, 0.03, 0.01).unwrap();
    for sampler in [Sampler::Iid, Sampler::OwenSobol] {
        let settings = settings(13, 9, 10, sampler);
        let baseline_policy = execution(13, 9, 1, 0x3000 + sampler as u64);
        let baseline = with_cx(|cx| {
            render_adaptive_with_execution(&scene, cx, &settings, policy, &baseline_policy)
        })
        .expect("serial-layout adaptive baseline");

        for (ordinal, (tile_width, tile_height, workers)) in
            [(4, 3, 2), (5, 2, 4), (7, 5, 8)].into_iter().enumerate()
        {
            let execution = execution(tile_width, tile_height, workers, 0x3100 + ordinal as u64)
                .with_quantum_weights((1..=workers as u32).collect())
                .expect("valid deterministic skew weights");
            let output = with_cx(|cx| {
                render_adaptive_with_execution(&scene, cx, &settings, policy, &execution)
            })
            .expect("parallel adaptive render");
            assert_adaptive_bits_eq(
                &baseline.film,
                &output.film,
                &format!("sampler={sampler:?} workers={workers} tile={tile_width}x{tile_height}"),
            );
        }

        let parked_policy = execution(5, 3, 4, 0x3200 + sampler as u64);
        let pool = RenderWorkerPool::new(&parked_policy, ExecMode::Deterministic, SEED);
        let parked = pool
            .with_parked_crew_local(|renderer| {
                with_cx(|cx| {
                    renderer.render_adaptive(&scene, cx, &settings, policy, &parked_policy)
                })
            })
            .expect("parked adaptive render");
        assert_adaptive_bits_eq(
            &baseline.film,
            &parked.film,
            &format!("sampler={sampler:?} parked crew"),
        );
    }
}

#[test]
fn g3_full_ceiling_preserves_the_uniform_raw_sum_for_both_samplers() {
    let scene = environment_scene();
    for sampler in [Sampler::Iid, Sampler::OwenSobol] {
        let settings = settings(12, 7, 6, sampler);
        let uniform = with_cx(|cx| render(&scene, cx, &settings)).expect("uniform oracle");
        let adaptive = with_cx(|cx| {
            render_adaptive_with_execution(
                &scene,
                cx,
                &settings,
                strict_policy(settings.spp),
                &execution(5, 3, 4, 0x4000 + sampler as u64),
            )
        })
        .expect("full-ceiling adaptive oracle");
        assert_uniform_raw_sum_bits_eq(
            &uniform,
            &adaptive.film,
            &format!("sampler={sampler:?} uniform raw sum"),
        );
        assert!(
            adaptive
                .film
                .sample_counts()
                .iter()
                .all(|samples| *samples == settings.spp)
        );
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
        "adaptive-battery-cancelling-sphere"
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
fn g4_owned_adaptive_row_prefix_resumes_to_uninterrupted_bits() {
    let settings = settings(35, 19, 3, Sampler::Iid);
    let policy = strict_policy(2);
    let execution = execution(8, 5, 4, 0x5001);
    let gate = Arc::new(CancelGate::new_clock_free());
    let evaluations = Arc::new(AtomicUsize::new(0));
    let resumable_scene = with_chart(
        environment_scene(),
        Box::new(CancellingSphere {
            evaluations: Arc::clone(&evaluations),
            cancel_at: Some(2_000),
            gate: Some(Arc::clone(&gate)),
        }),
    );
    let pending = with_cx(|cx| {
        PendingAdaptiveRender::begin_static(
            &resumable_scene,
            cx,
            settings,
            policy,
            execution.clone(),
        )
    })
    .expect("admit owned adaptive render");
    let suspended = with_gate_cx(&gate, |cx| pending.resume(cx))
        .expect_err("injected cancellation must suspend the adaptive job");
    assert!(matches!(
        suspended.cause(),
        RenderExecutionError::Tracer(TracerError::Cancelled)
    ));
    let partial = suspended.progress();
    assert!(
        partial.committed_tile_rows > 0 && partial.committed_tile_rows < partial.total_tile_rows,
        "cancellation must retain a strict complete-row prefix: {partial:?} evaluations={}",
        evaluations.load(Ordering::SeqCst)
    );
    assert_eq!(suspended.attempt_report().attempt_index, 1);
    assert!(suspended.attempt_report().memory.used_bytes > 0);

    let retry_pool = RenderWorkerPool::new(&execution, ExecMode::Deterministic, SEED);
    let resumed = retry_pool
        .with_parked_crew_local(|renderer| {
            with_cx(|cx| suspended.into_pending().resume_on_parked(renderer, cx))
        })
        .expect("fresh cancellation authority must finish exact adaptive resume");
    let reference_scene = with_chart(
        environment_scene(),
        Box::new(CancellingSphere {
            evaluations: Arc::new(AtomicUsize::new(0)),
            cancel_at: None,
            gate: None,
        }),
    );
    let reference = with_cx(|cx| {
        render_adaptive_with_execution(&reference_scene, cx, &settings, policy, &execution)
    })
    .expect("uninterrupted adaptive reference");
    assert_adaptive_bits_eq(
        &reference.film,
        &resumed.film,
        "cancel/resume adaptive exactness",
    );
    assert_eq!(resumed.report.attempt_index, 2);
    assert_eq!(resumed.report.memory.used_bytes, 0);
}
