//! Battery for the target-inaccessible ensemble executor
//! (bead frankensim-jmh21.2, core slice): coordinate-derived determinism,
//! complete accounting as a projection, rung admission, drain-on-cancel,
//! reserved-marker refusal, and the executor-to-bundle join.

use fs_blake3::ContentHash;
use fs_evidence::prediction_bundle::{
    AccessPolicy, ModelRungPolicy, OutputArtifactRef, OutputFamily, PredictionExecutionInput,
    PredictionOutputBundle, RandomStreamDesign, SampleAccounting,
};
use fs_evidence::vv::{ApplicabilityPolicy, ArtifactId, ArtifactKind, ArtifactRef};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_session::prediction_executor::{
    ExecutorRefusal, RunDisposition, SampleOutcome, execute_ensemble, sample_seeds,
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
                &[coordinates.sample_index.to_le_bytes().as_slice(), &seed.to_le_bytes()].concat(),
            )],
        }
    } else {
        SampleOutcome::Refused {
            rule: "test-parity-refusal".to_string(),
        }
    }
}

#[test]
fn seeds_are_pure_functions_of_logical_coordinates() {
    let input = admitted_input();
    // Recomputing in any order gives identical seeds; nothing about
    // execution context (workers, wall clock) participates.
    let forward: Vec<_> = (0..8).map(|index| sample_seeds(&input, index)).collect();
    let backward: Vec<_> = (0..8).rev().map(|index| sample_seeds(&input, index)).collect();
    for (index, seeds) in forward.iter().enumerate() {
        assert_eq!(seeds, &backward[7 - index]);
    }
    // Distinct samples get distinct seeds; distinct streams differ too.
    assert_ne!(forward[0].stream("sample-draw"), forward[1].stream("sample-draw"));
    assert_ne!(forward[0].stream("sample-draw"), forward[0].stream("jitter"));
    // Undeclared streams do not exist.
    assert_eq!(forward[0].stream("undeclared"), None);
}

#[test]
fn replay_is_bit_identical_and_accounting_is_a_projection() {
    let input = admitted_input();
    let first = with_cx(|cx, _| {
        execute_ensemble(cx, &input, "reduced-order", 64, parity_model).expect("runs")
    });
    let replay = with_cx(|cx, _| {
        execute_ensemble(cx, &input, "reduced-order", 64, parity_model).expect("runs")
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
        execute_ensemble(cx, &input, "reduced-order", 10, |coordinates, _| {
            if coordinates.sample_index == 3 {
                SampleOutcome::Failed {
                    rule: "test-numerical-failure".to_string(),
                }
            } else {
                SampleOutcome::Succeeded {
                    artifact_hashes: vec![fs_blake3::hash_bytes(b"artifact")],
                }
            }
        })
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
        execute_ensemble(cx, &input, "full-fem", 4, parity_model).expect_err("must refuse")
    });
    assert_eq!(error.rule, "prediction-executor-rung-not-admitted");
}

#[test]
fn ensemble_bounds_refuse_at_zero_and_cap_plus_one() {
    use fs_session::prediction_executor::MAX_ENSEMBLE_SAMPLES;
    let input = admitted_input();
    let zero = with_cx(|cx, _| execute_ensemble(cx, &input, "reduced-order", 0, parity_model));
    assert_eq!(zero.expect_err("zero refuses").rule, "prediction-executor-ensemble-bounds");
    let over = with_cx(|cx, _| {
        execute_ensemble(cx, &input, "reduced-order", MAX_ENSEMBLE_SAMPLES + 1, parity_model)
    });
    assert_eq!(over.expect_err("cap+1 refuses").rule, "prediction-executor-ensemble-bounds");
}

#[test]
fn cancellation_drains_marks_and_never_publishes_denominators() {
    let input = admitted_input();
    let run = with_cx(|cx, gate| {
        let mut executed = 0u64;
        execute_ensemble(cx, &input, "reduced-order", 100, |coordinates, seeds| {
            executed += 1;
            if executed == 5 {
                gate.cancel();
            }
            parity_model(coordinates, seeds)
        })
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
        execute_ensemble(cx, &input, "reduced-order", 2, |_, _| SampleOutcome::Cancelled)
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
