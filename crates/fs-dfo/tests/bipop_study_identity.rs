//! G3/G5 coverage for the production BIPOP full-study identity (7tv.23.4).
//!
//! The identity composes the exact root, callback content, retained trace
//! receipt, compatibility projections, and ordered restart ledger. These tests
//! exercise only public same-ISA replay/reference admission; private stale and
//! correctly resealed payload mutations live beside the sealed report fields.

#![deny(unsafe_code)]

use fs_dfo::{
    BIPOP_STUDY_IDENTITY_DOMAIN, BIPOP_STUDY_IDENTITY_SCHEMA_VERSION, BipopReport,
    BipopStudyAdmissionError, try_bipop_cmaes,
};

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
        "frankensim.fs-dfo.bipop-full-study.v1"
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
