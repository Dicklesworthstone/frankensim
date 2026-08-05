//! Focused G4/G5 battery for row-atomic voluntary yields and canonical
//! crash-recovery checkpoints. Failures report exact progress or binary64 bits.

#![cfg(feature = "tracer")]

use asupersync::types::Budget;
use fs_blake3::ContentHash;
use fs_exec::{CancelGate, Cx, ExecMode, RunId, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::lighting::EnvironmentMap;
use fs_render::tracer::{
    AdaptiveFilm, AdaptiveSamplingConfig, Camera, DirectStrategy, Film, PendingAdaptiveRender,
    PendingRender, RenderCheckpointBinding, RenderCheckpointError, RenderCheckpointKind,
    RenderCheckpointWriteError, RenderExecutionConfig, Sampler, Scene, Settings,
    render_adaptive_with_execution, render_with_execution,
};
use std::convert::Infallible;
use std::num::NonZeroU32;

const SEED: u64 = 0x6368_6563_6b70_7431;
const MAX_CHECKPOINT_BYTES: u64 = 8 << 20;

fn with_gate_budget_cx<R>(
    gate: &CancelGate,
    budget: Budget,
    operation: impl FnOnce(&Cx<'_>) -> R,
) -> R {
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: SEED,
                kernel_id: 0x4348_4b50,
                tile: 0,
                iteration: 0,
            },
            budget,
            ExecMode::Deterministic,
        );
        operation(&cx)
    })
}

fn with_gate_cx<R>(gate: &CancelGate, operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_gate_budget_cx(gate, Budget::INFINITE, operation)
}

fn with_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_gate_cx(&CancelGate::new_clock_free(), operation)
}

fn with_budget_cx<R>(budget: Budget, operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_gate_budget_cx(&CancelGate::new_clock_free(), budget, operation)
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
            .expect("admit deterministic environment"),
        ),
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 0.0),
            forward: Vec3::new(1.0, 0.0, 0.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.7,
        },
    }
}

fn settings(sampler: Sampler) -> Settings {
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

fn execution() -> RenderExecutionConfig {
    RenderExecutionConfig::try_new(3, 3, 2, 32 << 20, RunId(0x4350_0001))
        .expect("admit deterministic checkpoint execution")
}

fn adaptive_policy() -> AdaptiveSamplingConfig {
    AdaptiveSamplingConfig::try_new(2, 1, 0.0, 0.0, 0.0)
        .expect("admit deterministic adaptive policy")
}

fn identity(byte: u8) -> ContentHash {
    ContentHash([byte; 32])
}

fn binding_with_job(job: ContentHash) -> RenderCheckpointBinding {
    RenderCheckpointBinding::try_new(
        identity(1),
        identity(2),
        identity(3),
        identity(4),
        job,
        identity(6),
        identity(7),
        0,
        None,
    )
    .expect("admit complete root checkpoint binding")
}

fn write_uniform(
    pending: &PendingRender<'_>,
    binding: RenderCheckpointBinding,
) -> (Vec<u8>, fs_render::tracer::RenderCheckpointReceipt) {
    let mut bytes = Vec::new();
    let receipt = with_cx(|cx| {
        pending.write_checkpoint(
            binding,
            MAX_CHECKPOINT_BYTES,
            cx,
            |chunk| -> Result<(), Infallible> {
                bytes.extend_from_slice(chunk);
                Ok(())
            },
        )
    })
    .expect("stream uniform checkpoint");
    assert_eq!(receipt.byte_len(), bytes.len() as u64);
    (bytes, receipt)
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
                "{context}: pixel={pixel} channel={channel}"
            );
        }
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
    assert_eq!(left.time_mode(), right.time_mode(), "{context}");
    assert_eq!(left.sample_counts(), right.sample_counts(), "{context}");
    assert_eq!(left.decisions(), right.decisions(), "{context}");
    for pixel in 0..left.xyz_sums().len() {
        for channel in 0..3 {
            assert_eq!(
                left.xyz_sums()[pixel][channel].to_bits(),
                right.xyz_sums()[pixel][channel].to_bits(),
                "{context}: sum pixel={pixel} channel={channel}"
            );
            assert_eq!(
                left.running_means_xyz()[pixel][channel].to_bits(),
                right.running_means_xyz()[pixel][channel].to_bits(),
                "{context}: mean pixel={pixel} channel={channel}"
            );
            assert_eq!(
                left.m2_xyz()[pixel][channel].to_bits(),
                right.m2_xyz()[pixel][channel].to_bits(),
                "{context}: m2 pixel={pixel} channel={channel}"
            );
        }
    }
}

#[test]
fn g4_uniform_safe_point_checkpoint_restore_finishes_without_double_samples() {
    let settings = settings(Sampler::Iid);
    let execution = execution();
    let source_scene = scene();
    let pending =
        with_cx(|cx| PendingRender::begin_static(&source_scene, cx, settings, execution.clone()))
            .expect("begin uniform pending render");
    let yielded =
        with_cx(|cx| pending.advance_to_safe_point(cx, NonZeroU32::new(1).expect("nonzero quota")))
            .expect("intentional uniform checkpoint yield");
    let partial = yielded.progress();
    assert_eq!(partial.attempts, 1);
    assert_eq!(partial.committed_tile_rows, partial.total_tiles);
    assert!(
        partial.committed_tile_rows < partial.total_tile_rows,
        "{partial:?}"
    );
    assert_eq!(yielded.attempt_report().attempt_index, 1);

    let pending = yielded.into_pending();
    let binding = binding_with_job(pending.checkpoint_job_identity());
    let (bytes, written) = write_uniform(&pending, binding);
    assert_eq!(written.kind(), RenderCheckpointKind::Uniform);
    assert_eq!(written.progress(), partial);

    let restored_scene = scene();
    let (restored, read) = with_cx(|cx| {
        let fresh = PendingRender::begin_static(&restored_scene, cx, settings, execution.clone())
            .expect("re-admit uniform restore target");
        fresh.restore_checkpoint(binding, &bytes, MAX_CHECKPOINT_BYTES, cx)
    })
    .expect("restore uniform checkpoint");
    assert_eq!(read, written);
    assert_eq!(restored.progress(), partial);
    let resumed = with_cx(|cx| restored.resume(cx)).expect("finish restored uniform render");
    assert_eq!(resumed.report.attempt_index, 2);

    let reference_scene = scene();
    let reference =
        with_cx(|cx| render_with_execution(&reference_scene, cx, &settings, &execution))
            .expect("uninterrupted uniform reference");
    assert_film_bits_eq(
        &resumed.film,
        &reference.film,
        "uniform safe-point/checkpoint exactness",
    );
}

#[test]
fn g4_adaptive_safe_point_checkpoint_restore_finishes_with_exact_aovs() {
    let settings = settings(Sampler::OwenSobol);
    let policy = adaptive_policy();
    let execution = execution();
    let source_scene = scene();
    let pending = with_cx(|cx| {
        PendingAdaptiveRender::begin_static(&source_scene, cx, settings, policy, execution.clone())
    })
    .expect("begin adaptive pending render");
    let yielded =
        with_cx(|cx| pending.advance_to_safe_point(cx, NonZeroU32::new(1).expect("nonzero quota")))
            .expect("intentional adaptive checkpoint yield");
    let partial = yielded.progress();
    assert_eq!(partial.attempts, 1);
    assert_eq!(partial.committed_tile_rows, partial.total_tiles);
    assert!(
        partial.committed_tile_rows < partial.total_tile_rows,
        "{partial:?}"
    );
    assert_eq!(yielded.attempt_report().attempt_index, 1);

    let pending = yielded.into_pending();
    let binding = binding_with_job(pending.checkpoint_job_identity());
    let mut bytes = Vec::new();
    let written = with_cx(|cx| {
        pending.write_checkpoint(
            binding,
            MAX_CHECKPOINT_BYTES,
            cx,
            |chunk| -> Result<(), Infallible> {
                bytes.extend_from_slice(chunk);
                Ok(())
            },
        )
    })
    .expect("stream adaptive checkpoint");
    assert_eq!(written.kind(), RenderCheckpointKind::Adaptive);
    assert_eq!(written.progress(), partial);

    let restored_scene = scene();
    let (restored, read) = with_cx(|cx| {
        let fresh = PendingAdaptiveRender::begin_static(
            &restored_scene,
            cx,
            settings,
            policy,
            execution.clone(),
        )
        .expect("re-admit adaptive restore target");
        fresh.restore_checkpoint(binding, &bytes, MAX_CHECKPOINT_BYTES, cx)
    })
    .expect("restore adaptive checkpoint");
    assert_eq!(read, written);
    let completed = with_cx(|cx| {
        restored.advance_to_safe_point(
            cx,
            NonZeroU32::new(u32::MAX).expect("nonzero completion quota"),
        )
    })
    .expect("finish restored adaptive work at an opaque safe point");
    assert_eq!(completed.attempt_report().attempt_index, 2);
    assert_eq!(
        completed.progress().completed_tiles,
        completed.progress().total_tiles
    );
    assert_eq!(
        completed.attempt_report().executor.kernel,
        "fs-render/pending-adaptive-spectral-film-tile-v1"
    );

    let completed_again = with_cx(|cx| {
        completed
            .into_pending()
            .advance_to_safe_point(cx, NonZeroU32::new(1).expect("nonzero completed-job quota"))
    })
    .expect("yield an already-complete adaptive job without retracing");
    assert_eq!(completed_again.attempt_report().attempt_index, 3);
    assert_eq!(completed_again.attempt_report().workers, 0);
    assert_eq!(completed_again.attempt_report().executor.total, 0);
    assert_eq!(
        completed_again.attempt_report().executor.kernel,
        "fs-render/pending-adaptive-spectral-film-tile-v1"
    );
    let resumed = with_cx(|cx| completed_again.into_pending().resume(cx))
        .expect("publish completed restored adaptive render without retracing");
    assert_eq!(resumed.report.attempt_index, 4);
    assert_eq!(resumed.report.workers, 0);
    assert_eq!(
        resumed.report.executor.kernel,
        "fs-render/pending-adaptive-spectral-film-tile-v1"
    );

    let reference_scene = scene();
    let reference = with_cx(|cx| {
        render_adaptive_with_execution(&reference_scene, cx, &settings, policy, &execution)
    })
    .expect("uninterrupted adaptive reference");
    assert_adaptive_bits_eq(
        &resumed.film,
        &reference.film,
        "adaptive safe-point/checkpoint exactness",
    );
}

#[test]
fn g0_checkpoint_refuses_execution_budget_substitution_for_uniform_and_adaptive_jobs() {
    let settings = Settings {
        width: 2,
        height: 2,
        ..settings(Sampler::Iid)
    };
    let execution = execution();
    let source_scene = scene();

    let uniform =
        with_cx(|cx| PendingRender::begin_static(&source_scene, cx, settings, execution.clone()))
            .expect("admit infinite-budget uniform checkpoint source");
    let uniform_binding = binding_with_job(uniform.checkpoint_job_identity());
    let (uniform_bytes, _) = write_uniform(&uniform, uniform_binding);
    let restored_scene = scene();
    let uniform_refused = with_budget_cx(Budget::new().with_cost_quota(65_536), |cx| {
        let fresh = PendingRender::begin_static(&restored_scene, cx, settings, execution.clone())
            .expect("admit otherwise-identical finite-budget uniform restore target");
        fresh.restore_checkpoint(uniform_binding, &uniform_bytes, MAX_CHECKPOINT_BYTES, cx)
    });
    assert!(
        matches!(
            uniform_refused,
            Err(RenderCheckpointError::JobMismatch {
                field: "execution_budget"
            })
        ),
        "uniform restore accepted substituted finite execution budget: {uniform_refused:?}"
    );

    let policy = adaptive_policy();
    let adaptive = with_cx(|cx| {
        PendingAdaptiveRender::begin_static(&source_scene, cx, settings, policy, execution.clone())
    })
    .expect("admit infinite-budget adaptive checkpoint source");
    let adaptive_binding = binding_with_job(adaptive.checkpoint_job_identity());
    let mut adaptive_bytes = Vec::new();
    with_cx(|cx| {
        adaptive.write_checkpoint(
            adaptive_binding,
            MAX_CHECKPOINT_BYTES,
            cx,
            |chunk| -> Result<(), Infallible> {
                adaptive_bytes.extend_from_slice(chunk);
                Ok(())
            },
        )
    })
    .expect("stream infinite-budget adaptive checkpoint");
    let restored_scene = scene();
    let adaptive_refused = with_budget_cx(Budget::new().with_cost_quota(65_536), |cx| {
        let fresh = PendingAdaptiveRender::begin_static(
            &restored_scene,
            cx,
            settings,
            policy,
            execution.clone(),
        )
        .expect("admit otherwise-identical finite-budget adaptive restore target");
        fresh.restore_checkpoint(adaptive_binding, &adaptive_bytes, MAX_CHECKPOINT_BYTES, cx)
    });
    assert!(
        matches!(
            adaptive_refused,
            Err(RenderCheckpointError::JobMismatch {
                field: "execution_budget"
            })
        ),
        "adaptive restore accepted substituted finite execution budget: {adaptive_refused:?}"
    );
}

#[test]
fn g0_checkpoint_refuses_every_truncation_corruption_wrong_job_and_short_budget() {
    let settings = Settings {
        width: 2,
        height: 2,
        ..settings(Sampler::Iid)
    };
    let execution = execution();
    let source_scene = scene();
    let pending =
        with_cx(|cx| PendingRender::begin_static(&source_scene, cx, settings, execution.clone()))
            .expect("begin small checkpoint fixture");
    let binding = binding_with_job(pending.checkpoint_job_identity());
    let (bytes, receipt) = write_uniform(&pending, binding);

    let mut sink_called = false;
    let short_write = with_cx(|cx| {
        pending.write_checkpoint(
            binding,
            receipt.byte_len() - 1,
            cx,
            |_chunk| -> Result<(), Infallible> {
                sink_called = true;
                Ok(())
            },
        )
    })
    .expect_err("one-short encoder budget must refuse before emission");
    assert!(!sink_called);
    assert!(matches!(
        short_write,
        RenderCheckpointWriteError::Checkpoint(RenderCheckpointError::ByteLimitExceeded { .. })
    ));

    for prefix in 0..bytes.len() {
        let refused = with_cx(|cx| {
            let fresh = PendingRender::begin_static(&source_scene, cx, settings, execution.clone())
                .expect("re-admit truncation target");
            fresh.restore_checkpoint(binding, &bytes[..prefix], MAX_CHECKPOINT_BYTES, cx)
        });
        assert!(
            refused.is_err(),
            "accepted truncated prefix {prefix}/{}",
            bytes.len()
        );
    }

    for offset in [16, bytes.len() - 1] {
        let mut corrupt = bytes.clone();
        corrupt[offset] ^= 0x01;
        let refused = with_cx(|cx| {
            let fresh = PendingRender::begin_static(&source_scene, cx, settings, execution.clone())
                .expect("re-admit corruption target");
            fresh.restore_checkpoint(binding, &corrupt, MAX_CHECKPOINT_BYTES, cx)
        });
        assert!(
            matches!(refused, Err(RenderCheckpointError::IntegrityMismatch)),
            "corruption at byte {offset} produced {refused:?}"
        );
    }

    let wrong_job = with_cx(|cx| {
        let fresh = PendingRender::begin_static(&source_scene, cx, settings, execution.clone())
            .expect("re-admit wrong-job target");
        fresh.restore_checkpoint(
            binding_with_job(identity(99)),
            &bytes,
            MAX_CHECKPOINT_BYTES,
            cx,
        )
    });
    assert!(matches!(
        wrong_job,
        Err(RenderCheckpointError::BindingMismatch {
            field: "render_job_identity"
        })
    ));

    let short_read = with_cx(|cx| {
        let fresh = PendingRender::begin_static(&source_scene, cx, settings, execution.clone())
            .expect("re-admit short-read target");
        fresh.restore_checkpoint(binding, &bytes, bytes.len() as u64 - 1, cx)
    });
    assert!(matches!(
        short_read,
        Err(RenderCheckpointError::ByteLimitExceeded { .. })
    ));
}

#[test]
fn g4_checkpoint_observes_precancel_and_cancel_requested_by_final_seal_sink() {
    let settings = Settings {
        width: 2,
        height: 2,
        ..settings(Sampler::Iid)
    };
    let execution = execution();
    let source_scene = scene();
    let pending = with_cx(|cx| PendingRender::begin_static(&source_scene, cx, settings, execution))
        .expect("begin cancellation checkpoint fixture");
    let binding = binding_with_job(pending.checkpoint_job_identity());

    let precancelled = CancelGate::new_clock_free();
    precancelled.request();
    let mut pre_bytes = Vec::new();
    let pre = with_gate_cx(&precancelled, |cx| {
        pending.write_checkpoint(
            binding,
            MAX_CHECKPOINT_BYTES,
            cx,
            |chunk| -> Result<(), Infallible> {
                pre_bytes.extend_from_slice(chunk);
                Ok(())
            },
        )
    });
    assert!(matches!(
        pre,
        Err(RenderCheckpointWriteError::Checkpoint(
            RenderCheckpointError::Cancelled
        ))
    ));
    assert!(pre_bytes.is_empty());

    let gate = CancelGate::new_clock_free();
    let mut sealed_bytes = Vec::new();
    let after_seal = with_gate_cx(&gate, |cx| {
        pending.write_checkpoint(
            binding,
            MAX_CHECKPOINT_BYTES,
            cx,
            |chunk| -> Result<(), Infallible> {
                sealed_bytes.extend_from_slice(chunk);
                if chunk.starts_with(b"FSRSEAL1") {
                    gate.request();
                }
                Ok(())
            },
        )
    });
    assert!(matches!(
        after_seal,
        Err(RenderCheckpointWriteError::Checkpoint(
            RenderCheckpointError::Cancelled
        ))
    ));
    assert!(sealed_bytes.ends_with(&(sealed_bytes.len() as u64).to_le_bytes()));
}
