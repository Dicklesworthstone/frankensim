//! Battery for the target-inaccessible ensemble executor
//! (bead frankensim-jmh21.2): coordinate-derived determinism, complete
//! accounting as a projection, rung admission, drain-on-cancel,
//! reserved-marker refusal, the executor-to-bundle join, explicit
//! capability admission bound into checkpoint lineage and run logs,
//! honest fork lineage, the adaptive-diagnostics-only driver, and
//! worker-count-invariant replay.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_evidence::prediction_bundle::{
    AccessPolicy, ModelRungPolicy, OutputArtifactRef, OutputFamily, PredictionExecutionInput,
    PredictionOutputBundle, RandomStreamDesign, SampleAccounting,
};
use fs_evidence::vv::{ApplicabilityPolicy, ArtifactId, ArtifactKind, ArtifactRef};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_session::prediction_executor::{
    EnsembleCheckpoint, ExecutionCapabilities, ExecutorRefusal, RunDisposition, SampleOutcome,
    execute_ensemble, execute_ensemble_adaptive, fork_ensemble, resume_ensemble, sample_seeds,
};

fn reference(kind: ArtifactKind, id: &str, salt: u8) -> ArtifactRef {
    ArtifactRef::new(
        kind,
        ArtifactId::try_new(id).expect("valid id"),
        fs_blake3::hash_bytes(&[salt]),
    )
}

fn admitted_input() -> PredictionExecutionInput {
    PredictionExecutionInput::try_new(
        vec![("head_sha".to_string(), "deadbeef".to_string())],
        reference(ArtifactKind::ContextOfUse, "cou-1", 1),
        reference(ArtifactKind::ValidationPlan, "plan-1", 2),
        reference(ArtifactKind::CalibrationSplit, "split-1", 3),
        vec![reference(ArtifactKind::ExperimentArtifact, "scenario-a", 4)],
        Vec::new(),
        vec![
            RandomStreamDesign {
                name: "sample-draw".to_string(),
                seed_domain: "org.frankensim.test.sample-draw".to_string(),
                seed: 7,
                substreams: 2,
            },
            RandomStreamDesign {
                name: "jitter".to_string(),
                seed_domain: "org.frankensim.test.jitter".to_string(),
                seed: 11,
                substreams: 1,
            },
        ],
        ModelRungPolicy {
            allowed_rungs: vec!["reduced-order".to_string()],
            applicability: ApplicabilityPolicy::Refuse,
        },
        vec!["junction-maximum".to_string()],
        "blind-prediction".to_string(),
        reference(ArtifactKind::ExperimentArtifact, "holdout-1", 8),
        AccessPolicy::ExecutorOnly,
    )
    .expect("admitted input")
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x0000_C0FF_EE00_0002,
                kernel_id: 72,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

/// A deterministic model whose per-sample result is a pure function of its
/// coordinates and seeds: succeeds on even seed parity, refuses on odd.
fn parity_model(
    coordinates: &fs_session::prediction_executor::SampleCoordinates,
    seeds: &fs_session::prediction_executor::SampleSeeds,
) -> SampleOutcome {
    let seed = seeds.stream("sample-draw").expect("declared stream");
    if seed % 2 == 0 {
        SampleOutcome::Succeeded {
            artifact_hashes: vec![fs_blake3::hash_bytes(
                &[
                    coordinates.sample_index.to_le_bytes().as_slice(),
                    &seed.to_le_bytes(),
                ]
                .concat(),
            )],
        }
    } else {
        SampleOutcome::Refused {
            rule: "test-parity-refusal".to_string(),
        }
    }
}

/// The default grant of the pre-existing battery cases: compute only.
fn compute_only() -> ExecutionCapabilities {
    ExecutionCapabilities::compute_only()
}

/// Cancel a run right after `cancel_after` executed samples and return
/// its checkpoint — the standard staging fixture for resume/fork tests.
fn cancelled_checkpoint(
    input: &PredictionExecutionInput,
    requested: u64,
    cancel_after: u64,
    capabilities: ExecutionCapabilities,
) -> EnsembleCheckpoint {
    with_cx(|cx, gate| {
        let mut executed = 0u64;
        execute_ensemble(
            cx,
            input,
            "reduced-order",
            requested,
            capabilities,
            |coordinates, seeds| {
                executed += 1;
                if executed == cancel_after {
                    gate.request();
                }
                parity_model(coordinates, seeds)
            },
        )
        .expect("finalizes")
    })
    .checkpoint()
    .expect("cancelled runs checkpoint")
}

#[test]
fn seeds_are_pure_functions_of_logical_coordinates() {
    let input = admitted_input();
    // Recomputing in any order gives identical seeds; nothing about
    // execution context (workers, wall clock) participates.
    let forward: Vec<_> = (0..8).map(|index| sample_seeds(&input, index)).collect();
    let backward: Vec<_> = (0..8)
        .rev()
        .map(|index| sample_seeds(&input, index))
        .collect();
    for (index, seeds) in forward.iter().enumerate() {
        assert_eq!(seeds, &backward[7 - index]);
    }
    // Distinct samples get distinct seeds; distinct streams differ too.
    assert_ne!(
        forward[0].stream("sample-draw"),
        forward[1].stream("sample-draw")
    );
    assert_ne!(
        forward[0].stream("sample-draw"),
        forward[0].stream("jitter")
    );
    // Undeclared streams do not exist.
    assert_eq!(forward[0].stream("undeclared"), None);
}

#[test]
fn replay_is_bit_identical_and_accounting_is_a_projection() {
    let input = admitted_input();
    let first = with_cx(|cx, _| {
        execute_ensemble(
            cx,
            &input,
            "reduced-order",
            64,
            compute_only(),
            parity_model,
        )
        .expect("runs")
    });
    let replay = with_cx(|cx, _| {
        execute_ensemble(
            cx,
            &input,
            "reduced-order",
            64,
            compute_only(),
            parity_model,
        )
        .expect("runs")
    });
    assert_eq!(first, replay, "replay must be bit-identical");
    assert_eq!(first.disposition(), RunDisposition::Completed);

    let accounting = first.accounting().expect("completed run accounts");
    assert_eq!(accounting.requested, 64);
    assert_eq!(
        accounting.succeeded + accounting.refused + accounting.failed,
        64,
        "the partition is total"
    );
    // The projection agrees with the retained outcomes, by definition of
    // being derived from them - and both classes are genuinely populated.
    assert!(accounting.succeeded > 0 && accounting.refused > 0);
}

#[test]
fn every_failure_is_retained_and_reaches_the_denominator() {
    let input = admitted_input();
    let run = with_cx(|cx, _| {
        execute_ensemble(
            cx,
            &input,
            "reduced-order",
            10,
            compute_only(),
            |coordinates, _| {
                if coordinates.sample_index == 3 {
                    SampleOutcome::Failed {
                        rule: "test-numerical-failure".to_string(),
                    }
                } else {
                    SampleOutcome::Succeeded {
                        artifact_hashes: vec![fs_blake3::hash_bytes(b"artifact")],
                    }
                }
            },
        )
        .expect("runs")
    });
    assert!(matches!(run.outcomes()[3], SampleOutcome::Failed { .. }));
    let accounting = run.accounting().expect("accounts");
    assert_eq!((accounting.succeeded, accounting.failed), (9, 1));
    // The failed denominator flows into an admissible output bundle: the
    // scorer sees exactly 9/10, never 9/9.
    let bundle = PredictionOutputBundle::try_new(
        run.input_root(),
        accounting,
        vec![OutputArtifactRef {
            family: OutputFamily::Trajectory,
            id: "traj".to_string(),
            hash: fs_blake3::hash_bytes(b"traj"),
        }],
        Vec::new(),
        1,
        false,
        None,
        reference(ArtifactKind::SolutionVerificationReceipt, "sv-1", 9),
        vec!["recount".to_string()],
    )
    .expect("bundle admits the exact denominators");
    assert_eq!(bundle.accounting().failed, 1);
    bundle
        .verify_against_input(admitted_input().identity().expect("identity"))
        .expect("joins the sealed input");
}

#[test]
fn disallowed_rung_refuses_instead_of_substituting() {
    let input = admitted_input();
    let error = with_cx(|cx, _| {
        execute_ensemble(cx, &input, "full-fem", 4, compute_only(), parity_model)
            .expect_err("must refuse")
    });
    assert_eq!(error.rule, "prediction-executor-rung-not-admitted");
}

#[test]
fn ensemble_bounds_refuse_at_zero_and_cap_plus_one() {
    use fs_session::prediction_executor::MAX_ENSEMBLE_SAMPLES;
    let input = admitted_input();
    let zero = with_cx(|cx, _| {
        execute_ensemble(cx, &input, "reduced-order", 0, compute_only(), parity_model)
    });
    assert_eq!(
        zero.expect_err("zero refuses").rule,
        "prediction-executor-ensemble-bounds"
    );
    let over = with_cx(|cx, _| {
        execute_ensemble(
            cx,
            &input,
            "reduced-order",
            MAX_ENSEMBLE_SAMPLES + 1,
            compute_only(),
            parity_model,
        )
    });
    assert_eq!(
        over.expect_err("cap+1 refuses").rule,
        "prediction-executor-ensemble-bounds"
    );
}

#[test]
fn cancellation_drains_marks_and_never_publishes_denominators() {
    let input = admitted_input();
    let run = with_cx(|cx, gate| {
        let mut executed = 0u64;
        execute_ensemble(
            cx,
            &input,
            "reduced-order",
            100,
            compute_only(),
            |coordinates, seeds| {
                executed += 1;
                if executed == 5 {
                    gate.request();
                }
                parity_model(coordinates, seeds)
            },
        )
        .expect("cancelled runs still finalize")
    });
    let RunDisposition::Cancelled { drained_from } = run.disposition() else {
        panic!("expected cancellation, got {:?}", run.disposition());
    };
    assert_eq!(drained_from, 5, "poll happens before each sample");
    assert_eq!(run.outcomes().len(), 100, "every sample is present");
    assert!(
        run.outcomes()[usize::try_from(drained_from).expect("fits")..]
            .iter()
            .all(|outcome| *outcome == SampleOutcome::Cancelled),
        "every unexecuted sample carries the drain marker"
    );
    let error: ExecutorRefusal = run.accounting().expect_err("no partial denominators");
    assert_eq!(error.rule, "prediction-executor-cancelled-unscoreable");
}

#[test]
fn the_reserved_drain_marker_is_refused_from_models() {
    let input = admitted_input();
    let error = with_cx(|cx, _| {
        execute_ensemble(cx, &input, "reduced-order", 2, compute_only(), |_, _| {
            SampleOutcome::Cancelled
        })
        .expect_err("reserved marker refuses")
    });
    assert_eq!(error.rule, "prediction-executor-reserved-outcome");
}

#[test]
fn seed_derivation_is_domain_separated() {
    let input = admitted_input();
    let seeds = sample_seeds(&input, 0);
    let stream = &input.random_streams()[1]; // "sample-draw" sorts second? canonical order is by name
    // Recompute by hand through the published domain and formula.
    let mut payload = Vec::new();
    payload.extend_from_slice(stream.seed_domain.as_bytes());
    payload.push(0);
    payload.extend_from_slice(&stream.seed.to_le_bytes());
    payload.extend_from_slice(&0u64.to_le_bytes());
    let digest = fs_blake3::hash_domain(
        fs_session::prediction_executor::SAMPLE_SEED_DOMAIN,
        &payload,
    );
    let expected = u64::from_le_bytes(digest.as_bytes()[..8].try_into().expect("8 bytes"));
    assert_eq!(seeds.stream(&stream.name), Some(expected));
    // A bare (non-domain) hash of the same payload must NOT reproduce it.
    let bare: ContentHash = fs_blake3::hash_bytes(&payload);
    let bare_seed = u64::from_le_bytes(bare.as_bytes()[..8].try_into().expect("8 bytes"));
    assert_ne!(seeds.stream(&stream.name), Some(bare_seed));
}

#[test]
fn resume_equals_uninterrupted_bit_for_bit() {
    let input = admitted_input();
    let uninterrupted = with_cx(|cx, _| {
        execute_ensemble(
            cx,
            &input,
            "reduced-order",
            40,
            compute_only(),
            parity_model,
        )
        .expect("runs")
    });
    // Cancel after 12 executed samples, checkpoint, resume.
    let cancelled = with_cx(|cx, gate| {
        let mut executed = 0u64;
        execute_ensemble(
            cx,
            &input,
            "reduced-order",
            40,
            compute_only(),
            |coordinates, seeds| {
                executed += 1;
                if executed == 12 {
                    gate.request();
                }
                parity_model(coordinates, seeds)
            },
        )
        .expect("finalizes")
    });
    let checkpoint = cancelled.checkpoint().expect("cancelled runs checkpoint");
    assert_eq!(checkpoint.executed_len(), 12);
    let resumed = with_cx(|cx, _| {
        resume_ensemble(cx, &input, &checkpoint, 40, compute_only(), parity_model).expect("resumes")
    });
    assert_eq!(
        resumed, uninterrupted,
        "resume must reproduce the uninterrupted run exactly"
    );
    let accounting = resumed.accounting().expect("completed");
    assert_eq!(accounting.requested, 40);
}

#[test]
fn checkpoints_bind_lineage_and_refuse_foreign_resume() {
    let input = admitted_input();
    let cancelled = with_cx(|cx, gate| {
        let mut executed = 0u64;
        execute_ensemble(
            cx,
            &input,
            "reduced-order",
            10,
            compute_only(),
            |coordinates, seeds| {
                executed += 1;
                if executed == 3 {
                    gate.request();
                }
                parity_model(coordinates, seeds)
            },
        )
        .expect("finalizes")
    });
    let checkpoint = cancelled.checkpoint().expect("checkpoints");

    // A completed run has nothing to resume.
    let completed = with_cx(|cx, _| {
        execute_ensemble(cx, &input, "reduced-order", 4, compute_only(), parity_model)
            .expect("runs")
    });
    assert_eq!(
        completed.checkpoint().expect_err("nothing to resume").rule,
        "prediction-executor-nothing-to-resume"
    );

    // Foreign input root refuses: build a second admitted input differing
    // in one stream seed.
    let mut streams = input.random_streams().to_vec();
    streams[1].seed ^= 1;
    let foreign = PredictionExecutionInput::try_new(
        vec![("head_sha".to_string(), "deadbeef".to_string())],
        reference(ArtifactKind::ContextOfUse, "cou-1", 1),
        reference(ArtifactKind::ValidationPlan, "plan-1", 2),
        reference(ArtifactKind::CalibrationSplit, "split-1", 3),
        vec![reference(ArtifactKind::ExperimentArtifact, "scenario-a", 4)],
        Vec::new(),
        streams,
        ModelRungPolicy {
            allowed_rungs: vec!["reduced-order".to_string()],
            applicability: ApplicabilityPolicy::Refuse,
        },
        vec!["junction-maximum".to_string()],
        "blind-prediction".to_string(),
        reference(ArtifactKind::ExperimentArtifact, "holdout-1", 8),
        AccessPolicy::ExecutorOnly,
    )
    .expect("admits");
    let error = with_cx(|cx, _| {
        resume_ensemble(cx, &foreign, &checkpoint, 10, compute_only(), parity_model)
            .expect_err("foreign root refuses")
    });
    assert_eq!(error.rule, "prediction-executor-foreign-checkpoint");

    // A prefix that is not strictly shorter than the request refuses.
    let error = with_cx(|cx, _| {
        resume_ensemble(cx, &input, &checkpoint, 3, compute_only(), parity_model)
            .expect_err("non-prefix refuses")
    });
    assert_eq!(error.rule, "prediction-executor-checkpoint-bounds");

    // The lineage identity moves with any retained outcome, the rung, or
    // the root - a tampered prefix is attestably different.
    let identity = checkpoint.identity();
    let retampered = with_cx(|cx, gate| {
        let mut executed = 0u64;
        execute_ensemble(
            cx,
            &input,
            "reduced-order",
            10,
            compute_only(),
            |coordinates, seeds| {
                executed += 1;
                if executed == 3 {
                    gate.request();
                }
                // Same cancellation point, different outcome content.
                match parity_model(coordinates, seeds) {
                    SampleOutcome::Refused { .. } => SampleOutcome::Refused {
                        rule: "different-rule".to_string(),
                    },
                    other => other,
                }
            },
        )
        .expect("finalizes")
    });
    let other = retampered.checkpoint().expect("checkpoints");
    if other.executed_len() == checkpoint.executed_len() && other != checkpoint {
        assert_ne!(
            other.identity(),
            identity,
            "outcome edits must move lineage"
        );
    }
    // Trivial self-consistency: identical checkpoints share identity.
    assert_eq!(checkpoint.identity(), identity);
}

#[test]
fn run_log_is_deterministic_bounded_and_content_addressed() {
    let input = admitted_input();
    let run = |_: ()| {
        with_cx(|cx, _| {
            execute_ensemble(
                cx,
                &input,
                "reduced-order",
                32,
                compute_only(),
                parity_model,
            )
            .expect("runs")
        })
    };
    let first = run(()).log();
    let replay = run(()).log();
    assert_eq!(
        first.canonical_bytes(),
        replay.canonical_bytes(),
        "replayed logs must be byte-identical"
    );
    assert_eq!(first.identity(), replay.identity());
    // Bounded: one class byte per sample plus deduplicated rules.
    assert!(
        first.canonical_bytes().len() < 32 + 2048,
        "log stays bounded"
    );
    // The reproduction command is repository-relative: no absolute paths.
    let reproduction = first.reproduction_command();
    assert!(reproduction.starts_with("cargo test -p fs-session"));
    assert!(!reproduction.contains("/Users") && !reproduction.contains("/home"));
    // Redaction by construction: the canonical bytes of two logs from the
    // same logical run never differ, so nothing environment-dependent
    // (wall clock, PID, host, worker) can be present.
}

#[test]
fn run_log_first_divergence_and_rule_counts_are_exact() {
    let input = admitted_input();
    let run = with_cx(|cx, _| {
        execute_ensemble(
            cx,
            &input,
            "reduced-order",
            12,
            compute_only(),
            |coordinates, _| match coordinates.sample_index {
                4 => SampleOutcome::Failed {
                    rule: "test-blowup".to_string(),
                },
                7 | 9 => SampleOutcome::Refused {
                    rule: "test-out-of-domain".to_string(),
                },
                _ => SampleOutcome::Succeeded {
                    artifact_hashes: vec![fs_blake3::hash_bytes(b"a")],
                },
            },
        )
        .expect("runs")
    });
    let log = run.log();
    assert_eq!(log.first_divergence(), Some(4));
    // Identity moves with outcome content: a run differing only in one
    // sample's rule has a different log identity.
    let other = with_cx(|cx, _| {
        execute_ensemble(
            cx,
            &input,
            "reduced-order",
            12,
            compute_only(),
            |coordinates, _| match coordinates.sample_index {
                4 => SampleOutcome::Failed {
                    rule: "test-DIFFERENT-blowup".to_string(),
                },
                7 | 9 => SampleOutcome::Refused {
                    rule: "test-out-of-domain".to_string(),
                },
                _ => SampleOutcome::Succeeded {
                    artifact_hashes: vec![fs_blake3::hash_bytes(b"a")],
                },
            },
        )
        .expect("runs")
    });
    assert_ne!(other.log().identity(), log.identity());
}

#[test]
fn capabilities_bind_lineage_and_resume_cannot_change_grant() {
    let input = admitted_input();
    let checkpoint = cancelled_checkpoint(&input, 10, 3, compute_only());
    assert_eq!(checkpoint.capabilities(), compute_only());

    // Resuming under a broader grant refuses...
    let escalated = with_cx(|cx, _| {
        resume_ensemble(
            cx,
            &input,
            &checkpoint,
            10,
            compute_only().granting_filesystem(),
            parity_model,
        )
    });
    assert_eq!(
        escalated.expect_err("escalation refuses").rule,
        "prediction-executor-capability-mismatch"
    );
    // ...and so does ANY other change of grant, including a different
    // same-width one: resume continues THE run, not a lookalike.
    let swapped = with_cx(|cx, _| {
        resume_ensemble(
            cx,
            &input,
            &checkpoint,
            10,
            compute_only().granting_network(),
            parity_model,
        )
    });
    assert_eq!(
        swapped.expect_err("any mismatch refuses").rule,
        "prediction-executor-capability-mismatch"
    );

    // The exact grant resumes to the uninterrupted run bit-for-bit.
    let uninterrupted = with_cx(|cx, _| {
        execute_ensemble(
            cx,
            &input,
            "reduced-order",
            10,
            compute_only(),
            parity_model,
        )
        .expect("runs")
    });
    let resumed = with_cx(|cx, _| {
        resume_ensemble(cx, &input, &checkpoint, 10, compute_only(), parity_model).expect("resumes")
    });
    assert_eq!(resumed, uninterrupted);

    // The lineage identity moves with the capability set itself.
    let privileged = cancelled_checkpoint(&input, 10, 3, compute_only().granting_filesystem());
    assert_ne!(
        privileged.identity(),
        checkpoint.identity(),
        "capability changes must move checkpoint lineage"
    );
}

#[test]
fn fork_records_lineage_and_cannot_pose_as_a_plain_run() {
    let input = admitted_input();
    let plain = with_cx(|cx, _| {
        execute_ensemble(
            cx,
            &input,
            "reduced-order",
            40,
            compute_only(),
            parity_model,
        )
        .expect("runs")
    });
    let checkpoint = cancelled_checkpoint(&input, 40, 12, compute_only());

    // A same-grant fork reproduces the exact outcome CONTENT of the plain
    // ensemble - determinism does not care how the prefix was staged -
    let forked = with_cx(|cx, _| {
        fork_ensemble(cx, &input, &checkpoint, 40, compute_only(), parity_model).expect("forks")
    });
    assert_eq!(forked.outcomes(), plain.outcomes());
    assert_eq!(forked.disposition(), RunDisposition::Completed);
    assert_eq!(forked.accounting().expect("completed").requested, 40);
    // ...but it can never POSE as that ensemble: the fork carries its
    // parent, and lineage participates in both run and log identity.
    assert_eq!(forked.fork_parent(), Some(checkpoint.identity()));
    assert_eq!(plain.fork_parent(), None);
    assert_ne!(forked, plain);
    assert_ne!(forked.log().identity(), plain.log().identity());

    // A fork under a BROADER grant is legitimate because it is recorded:
    // parent identity AND child grant both bind the child's log.
    let privileged = with_cx(|cx, _| {
        fork_ensemble(
            cx,
            &input,
            &checkpoint,
            40,
            compute_only().granting_filesystem(),
            parity_model,
        )
        .expect("forks")
    });
    assert_eq!(privileged.fork_parent(), Some(checkpoint.identity()));
    assert!(privileged.capabilities().admits_filesystem());
    assert!(!forked.capabilities().admits_filesystem());
    assert_ne!(
        privileged.log().identity(),
        forked.log().identity(),
        "the grant binds the forked run's log"
    );
}

#[test]
fn adaptive_history_exposes_every_retained_outcome_and_changes_nothing() {
    let input = admitted_input();
    let observed_lens = std::cell::RefCell::new(Vec::new());
    let saw_planted_failure = std::cell::Cell::new(false);
    // One shared decision function: the adaptive closure observes history
    // but MUST NOT change behavior, so the metamorphic comparison below is
    // against the SAME per-sample logic through the plain executor.
    fn decide(
        coordinates: &fs_session::prediction_executor::SampleCoordinates,
        seeds: &fs_session::prediction_executor::SampleSeeds,
    ) -> SampleOutcome {
        match coordinates.sample_index {
            3 => SampleOutcome::Failed {
                rule: "test-blowup".to_string(),
            },
            _ => parity_model(coordinates, seeds),
        }
    }
    let adaptive = with_cx(|cx, _| {
        execute_ensemble_adaptive(
            cx,
            &input,
            "reduced-order",
            12,
            compute_only(),
            |coordinates, seeds, history| {
                observed_lens.borrow_mut().push(history.len() as u64);
                // The retained prefix arrives VERBATIM: the failure planted
                // at index 3 stays visible to every later sample. There is
                // no API through which it could be dropped.
                if coordinates.sample_index > 3
                    && matches!(history[3], SampleOutcome::Failed { .. })
                {
                    saw_planted_failure.set(true);
                }
                decide(coordinates, seeds)
            },
        )
        .expect("runs")
    });
    assert_eq!(
        *observed_lens.borrow(),
        (0..12).collect::<Vec<u64>>(),
        "each sample sees exactly its own retained prefix"
    );
    assert!(
        saw_planted_failure.get(),
        "later samples observe the planted failure"
    );
    // Ignoring the history is EXACTLY the plain executor.
    let plain = with_cx(|cx, _| {
        execute_ensemble(
            cx,
            &input,
            "reduced-order",
            12,
            compute_only(),
            |coordinates, seeds| decide(coordinates, seeds),
        )
        .expect("runs")
    });
    assert_eq!(adaptive, plain);
}

#[test]
fn worker_partition_replays_bit_identically_for_every_chunking() {
    use std::cell::Cell;
    let input = admitted_input();
    let total = 24u64;
    let uninterrupted = with_cx(|cx, _| {
        execute_ensemble(
            cx,
            &input,
            "reduced-order",
            total,
            compute_only(),
            parity_model,
        )
        .expect("runs")
    });
    // Simulate W workers executing disjoint chunks sequentially: at each
    // chunk boundary the current stage cancels, checkpoints, and hands
    // the remainder to the next stage through the SAME coordinate-derived
    // seeds. Every partition must land on the identical run.
    for workers in [1u64, 2, 3, 4, 6, 8, 12] {
        let chunk = total / workers;
        let mut staged = with_cx(|cx, gate| {
            let next_cancel = Cell::new(chunk);
            execute_ensemble(
                cx,
                &input,
                "reduced-order",
                total,
                compute_only(),
                |coordinates, seeds| {
                    if coordinates.sample_index + 1 == next_cancel.get()
                        && next_cancel.get() < total
                    {
                        gate.request();
                    }
                    if coordinates.sample_index + 1 == next_cancel.get() {
                        next_cancel.set(next_cancel.get() + chunk);
                    }
                    parity_model(coordinates, seeds)
                },
            )
            .expect("stage finalizes")
        });
        while staged.disposition() != RunDisposition::Completed {
            let checkpoint = staged.checkpoint().expect("staged checkpoints");
            staged = with_cx(|cx, gate| {
                let next_cancel = Cell::new(checkpoint.executed_len() + chunk);
                resume_ensemble(
                    cx,
                    &input,
                    &checkpoint,
                    total,
                    compute_only(),
                    |coordinates, seeds| {
                        if coordinates.sample_index + 1 == next_cancel.get()
                            && next_cancel.get() < total
                        {
                            gate.request();
                        }
                        if coordinates.sample_index + 1 == next_cancel.get() {
                            next_cancel.set(next_cancel.get() + chunk);
                        }
                        parity_model(coordinates, seeds)
                    },
                )
                .expect("stage resumes")
            });
        }
        assert_eq!(
            staged, uninterrupted,
            "{workers}-worker partition replays bit-identically"
        );
    }
}

#[test]
fn run_log_v2_binds_capabilities_and_fork_parent() {
    let input = admitted_input();
    let plain = with_cx(|cx, _| {
        execute_ensemble(
            cx,
            &input,
            "reduced-order",
            16,
            compute_only(),
            parity_model,
        )
        .expect("runs")
    });
    let privileged = with_cx(|cx, _| {
        execute_ensemble(
            cx,
            &input,
            "reduced-order",
            16,
            compute_only().granting_filesystem(),
            parity_model,
        )
        .expect("runs")
    });
    // Identical outcomes, different declared environment: the LOG must
    // distinguish them even though nothing else about the bytes differs.
    assert_eq!(plain.outcomes(), privileged.outcomes());
    assert_ne!(plain.log().identity(), privileged.log().identity());
    assert!(
        fs_session::prediction_executor::RUN_LOG_SCHEMA.ends_with(".v2"),
        "the capability binding is a schema-level version bump"
    );
    // A fork is logged as a fork even when its shape matches a plain run.
    let checkpoint = cancelled_checkpoint(&input, 16, 8, compute_only());
    let fork_log_identity = with_cx(|cx, _| {
        fork_ensemble(cx, &input, &checkpoint, 16, compute_only(), parity_model)
            .expect("forks")
            .log()
            .identity()
    });
    assert_ne!(fork_log_identity, plain.log().identity());
}
