//! G3/G5 coverage for the production BIPOP full-study identity (7tv.23.4).
//!
//! The identity composes the exact root, callback content, retained trace
//! receipt, compatibility projections, and ordered restart ledger. These tests
//! exercise only public same-ISA replay/reference admission; private stale and
//! correctly resealed payload mutations live beside the sealed report fields.

#![deny(unsafe_code)]

use fs_dfo::{
    BIPOP_STUDY_IDENTITY_DOMAIN, BIPOP_STUDY_IDENTITY_SCHEMA_VERSION, BipopError,
    BipopReplayAdmissionError, BipopReport, BipopStudyAdmissionError, try_bipop_cmaes,
};
use std::cell::Cell;

const ROOT_SEED: u64 = 0xB1_90_00_04;

fn sphere(point: &[f64]) -> f64 {
    point.iter().map(|value| value * value).sum()
}

fn run_study(seed: u64) -> BipopReport {
    let mut objective = sphere;
    try_bipop_cmaes(&mut objective, &[2.0, -1.0], 0.75, 20, None, seed)
        .expect("fixed full-study identity fixture admits")
}

/// G5: same root and same callback semantics replay to the exact complete
/// production identity, and either replay admits against the retained reference.
#[test]
fn g5_same_seed_replays_and_admits_the_full_study_identity() {
    let first = run_study(ROOT_SEED);
    let replay = run_study(ROOT_SEED);
    first.validate_ledger().expect("first ledger validates");
    replay.validate_ledger().expect("replay ledger validates");
    first
        .validate_study_identity()
        .expect("first retained identity matches its payload");
    replay
        .validate_study_identity()
        .expect("replay retained identity matches its payload");

    let identity = first.study_identity();
    assert_eq!(identity, replay.study_identity());
    assert_eq!(
        identity.schema_version(),
        BIPOP_STUDY_IDENTITY_SCHEMA_VERSION
    );
    assert_eq!(identity.restarts(), first.records().len());
    assert_eq!(identity.evaluations(), first.total_evals);
    assert_ne!(identity.digest(), &[0u8; 32]);
    assert_eq!(
        BIPOP_STUDY_IDENTITY_DOMAIN,
        "frankensim.fs-dfo.bipop-full-study.v2"
    );
    first
        .admit_study_identity(identity)
        .expect("canonical payload admits against itself");
    replay
        .admit_study_identity(identity)
        .expect("same-seed replay admits against the canonical reference");
}

/// G3: a different causal seed remains internally valid but names a different
/// study, so reference admission must not misclassify it as a stale payload.
#[test]
fn g3_valid_different_seed_is_a_reference_identity_mismatch() {
    let reference = run_study(ROOT_SEED).study_identity();
    let other = run_study(ROOT_SEED + 1);
    other
        .validate_study_identity()
        .expect("different-seed payload is internally self-consistent");
    let found = other.study_identity();
    assert_ne!(found, reference);
    assert_eq!(
        other.admit_study_identity(reference),
        Err(BipopStudyAdmissionError::ReferenceIdentityMismatch {
            expected: reference,
            found,
        })
    );
}

/// G5: replay-backed admission invokes an independently supplied objective only
/// after cheap evidence admission succeeds, then requires the full replayed
/// production study identity rather than a summary-field match.
#[test]
fn g5_replay_backed_admission_accepts_the_same_semantic_study() {
    let report = run_study(ROOT_SEED);
    let identity = report.study_identity();
    let calls = Cell::new(0usize);
    let mut oracle = |point: &[f64]| {
        calls.set(calls.get() + 1);
        sphere(point)
    };

    report
        .admit_study_identity_with_replay(identity, &mut oracle)
        .expect("same root and objective replay to the complete retained study");
    assert_eq!(calls.get(), report.total_evals);
}

/// G3/G5: stale payloads and wrong external references are callback-free even
/// at the replay-backed surface. Classification is inherited exactly from the
/// cheap gate rather than being masked by a newly executed oracle.
#[test]
fn g3_replay_backed_admission_keeps_stale_and_reference_refusals_callback_free() {
    let canonical = run_study(ROOT_SEED);
    let reference = canonical.study_identity();

    let mut stale = canonical.clone();
    stale.best.f_best = f64::from_bits(stale.best.f_best.to_bits() ^ 1);
    let stale_calls = Cell::new(0usize);
    let mut stale_oracle = |point: &[f64]| {
        stale_calls.set(stale_calls.get() + 1);
        sphere(point)
    };
    let stale_error = stale
        .admit_study_identity_with_replay(reference, &mut stale_oracle)
        .expect_err("stale payload must refuse before replay");
    assert!(matches!(
        stale_error,
        BipopReplayAdmissionError::EvidenceAdmission {
            error: BipopStudyAdmissionError::PayloadIdentityMismatch { .. }
        }
    ));
    assert_eq!(stale_calls.get(), 0);

    let other = run_study(ROOT_SEED + 1);
    let other_calls = Cell::new(0usize);
    let mut other_oracle = |point: &[f64]| {
        other_calls.set(other_calls.get() + 1);
        sphere(point)
    };
    assert!(matches!(
        other.admit_study_identity_with_replay(reference, &mut other_oracle),
        Err(BipopReplayAdmissionError::EvidenceAdmission {
            error: BipopStudyAdmissionError::ReferenceIdentityMismatch { .. }
        })
    ));
    assert_eq!(other_calls.get(), 0);
}

/// G3: a different executable objective is a semantic mismatch after a complete
/// valid replay, not a stale-payload or external-reference classification.
#[test]
fn g3_replay_backed_admission_types_completed_semantic_mismatch() {
    let report = run_study(ROOT_SEED);
    let retained = report.study_identity();
    let calls = Cell::new(0usize);
    let mut different_oracle = |point: &[f64]| {
        calls.set(calls.get() + 1);
        sphere(point) + 1.0
    };

    let error = report
        .admit_study_identity_with_replay(retained, &mut different_oracle)
        .expect_err("different objective semantics must move complete replay evidence");
    match error {
        BipopReplayAdmissionError::SemanticMismatch {
            retained: found_retained,
            replayed,
        } => {
            assert_eq!(found_retained, retained);
            assert_ne!(replayed, retained);
        }
        other => panic!("expected completed semantic mismatch, found {other:?}"),
    }
    assert!(calls.get() > 0);
}

/// G4: replay execution refusal is distinct from a completed semantic mismatch.
/// The retained study stops safely at its initial target callback; the alternate
/// oracle misses that target and deterministically enters an overflowing first
/// candidate under the retained extreme root.
#[test]
fn g4_replay_backed_admission_types_execution_refusal() {
    const OVERFLOW_SEED: u64 = 0xB1_90_00_02;
    let mut at_target = |_point: &[f64]| 0.0;
    let report = try_bipop_cmaes(
        &mut at_target,
        &[f64::MAX],
        f64::MAX,
        5,
        Some(0.0),
        OVERFLOW_SEED,
    )
    .expect("initial target stops before the extreme root generates a candidate");
    let identity = report.study_identity();
    let calls = Cell::new(0usize);
    let mut misses_target = |_point: &[f64]| {
        calls.set(calls.get() + 1);
        1.0
    };

    assert!(matches!(
        report.admit_study_identity_with_replay(identity, &mut misses_target),
        Err(BipopReplayAdmissionError::ReplayExecution {
            error: BipopError::GeneratedCandidateNonFinite {
                restart: 0,
                generation: 1,
                candidate: 0,
                component: 0,
                ..
            }
        })
    ));
    assert_eq!(
        calls.get(),
        1,
        "only the finite initial replay point is observed"
    );
}
