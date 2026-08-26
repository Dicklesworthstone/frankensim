//! Adversarial Gauntlet and independent checker verification suite (bead `frankensim-mkfvu.3`).
//!
//! Tests:
//! - G0 algebraic laws & independent key reconstruction;
//! - G3 metamorphic attacks on field permutations and normalization;
//! - G4 tamper detection on declared vs computed artifact hashes (`InvalidEvidence`);
//! - G5 determinism verification and divergence refutation (`RefutedDivergence`);
//! - Nondeterministic mode handling (`ExplicitlyNondeterministic`).

use std::collections::BTreeMap;

use fs_ledger::{ContentHash, hash_bytes};
use fs_recompute::checker::{CheckerDisposition, IndependentChecker};
use fs_recompute::semantic_determinism::{
    ComputationKey, DeterminismClass, ExecutionPolicy, OutputObservation, ToleranceRole,
};

fn sample_code() -> ContentHash {
    hash_bytes(b"fs-solver-v0.2.0-kernel")
}

fn sample_input(id: u8) -> ContentHash {
    hash_bytes(&[id; 32])
}

#[test]
fn gauntlet_001_independent_key_hash_parity() {
    let code = sample_code();
    let mut params = BTreeMap::new();
    params.insert("tol".to_string(), "1e-5".to_string());
    params.insert("max_iter".to_string(), "100".to_string());

    let policy = ExecutionPolicy::try_new(
        DeterminismClass::ToleranceDependentDeterministic,
        ToleranceRole::StoppingCriterion,
        Some(1e-5),
        777,
        code,
        Some(100),
    )
    .unwrap();

    let key = ComputationKey::try_new(
        "iterative_solve",
        vec![sample_input(1), sample_input(2)],
        params.clone(),
        &policy,
    )
    .unwrap();

    let raw_params: Vec<(String, String)> = params.into_iter().collect();

    let independent_hash = IndependentChecker::compute_key_hash_independent(
        "iterative_solve",
        &[sample_input(1), sample_input(2)],
        &raw_params,
        DeterminismClass::ToleranceDependentDeterministic,
        ToleranceRole::StoppingCriterion,
        (1e-5f64).to_bits(),
        777,
        &code,
        Some(100),
    );

    assert_eq!(key.content_hash(), independent_hash);
}

#[test]
fn gauntlet_002_verified_policy_match_on_repeat() {
    let mut checker = IndependentChecker::new();
    let code = sample_code();
    let policy = ExecutionPolicy::exact_deterministic(code, 42);
    let key = ComputationKey::try_new("exact_op", vec![sample_input(1)], BTreeMap::new(), &policy)
        .unwrap();

    let artifact = b"exact-artifact-data";
    let art_hash = fs_recompute::artifact_content_hash(artifact);
    let obs = OutputObservation::try_new(art_hash, Some(0.0), Some(0.01), Some(1024)).unwrap();

    // First run -> VerifiedPolicyMatch
    let res1 = checker.check_observation(&key, &obs, artifact);
    assert_eq!(
        res1,
        CheckerDisposition::VerifiedPolicyMatch {
            computation_hash: key.content_hash(),
            artifact_hash: art_hash,
        }
    );

    // Repeat identical run -> VerifiedPolicyMatch
    let res2 = checker.check_observation(&key, &obs, artifact);
    assert_eq!(
        res2,
        CheckerDisposition::VerifiedPolicyMatch {
            computation_hash: key.content_hash(),
            artifact_hash: art_hash,
        }
    );
}

#[test]
fn gauntlet_003_divergence_refuted_with_exact_witness() {
    let mut checker = IndependentChecker::new();
    let code = sample_code();
    let policy = ExecutionPolicy::exact_deterministic(code, 42);
    let key = ComputationKey::try_new("exact_op", vec![sample_input(1)], BTreeMap::new(), &policy)
        .unwrap();

    let artifact1 = b"exact-artifact-data-1";
    let art_hash1 = fs_recompute::artifact_content_hash(artifact1);
    let obs1 = OutputObservation::try_new(art_hash1, Some(0.0), Some(0.01), Some(1024)).unwrap();
    checker.check_observation(&key, &obs1, artifact1);

    // Divergent run -> RefutedDivergence
    let artifact2 = b"exact-artifact-data-2";
    let art_hash2 = fs_recompute::artifact_content_hash(artifact2);
    let obs2 = OutputObservation::try_new(art_hash2, Some(0.0), Some(0.01), Some(1024)).unwrap();

    let res_div = checker.check_observation(&key, &obs2, artifact2);
    assert_eq!(
        res_div,
        CheckerDisposition::RefutedDivergence {
            computation_hash: key.content_hash(),
            expected: art_hash1,
            observed: art_hash2,
            first_divergent_field: "artifact_content_bytes",
        }
    );
    assert!(res_div.is_refuted());
}

#[test]
fn gauntlet_004_invalid_evidence_on_tampered_declared_hash() {
    let mut checker = IndependentChecker::new();
    let code = sample_code();
    let policy = ExecutionPolicy::exact_deterministic(code, 42);
    let key = ComputationKey::try_new("exact_op", vec![sample_input(1)], BTreeMap::new(), &policy)
        .unwrap();

    let artifact = b"actual-bytes";
    let fake_art_hash = hash_bytes(b"tampered-fake-hash");
    let obs = OutputObservation::try_new(fake_art_hash, None, None, None).unwrap();

    let res = checker.check_observation(&key, &obs, artifact);
    assert!(matches!(res, CheckerDisposition::InvalidEvidence { .. }));
}

#[test]
fn gauntlet_005_nondeterministic_mode_never_claims_divergence() {
    let mut checker = IndependentChecker::new();
    let code = sample_code();
    let policy = ExecutionPolicy::try_new(
        DeterminismClass::Nondeterministic,
        ToleranceRole::None,
        None,
        999,
        code,
        None,
    )
    .unwrap();
    let key = ComputationKey::try_new("heuristic", vec![], BTreeMap::new(), &policy).unwrap();

    let art1 = b"fast-bytes-1";
    let obs1 =
        OutputObservation::try_new(fs_recompute::artifact_content_hash(art1), None, None, None)
            .unwrap();
    let res1 = checker.check_observation(&key, &obs1, art1);
    assert_eq!(
        res1,
        CheckerDisposition::ExplicitlyNondeterministic {
            computation_hash: key.content_hash(),
        }
    );

    let art2 = b"fast-bytes-2";
    let obs2 =
        OutputObservation::try_new(fs_recompute::artifact_content_hash(art2), None, None, None)
            .unwrap();
    let res2 = checker.check_observation(&key, &obs2, art2);
    assert_eq!(
        res2,
        CheckerDisposition::ExplicitlyNondeterministic {
            computation_hash: key.content_hash(),
        }
    );
}

/// Regression (fresh-eyes review, 2026-08-26): the checker is the
/// independent-validation boundary; an observation carrying NaN evidence —
/// constructible via struct literal bypassing OutputObservation::try_new —
/// must receive InvalidEvidence and must never enter history.
#[test]
fn gauntlet_nan_observation_is_invalid_evidence() {
    use fs_recompute::NodeRecord;
    let code = sample_code();
    let policy = ExecutionPolicy::try_new(
        DeterminismClass::ToleranceDependentDeterministic,
        ToleranceRole::StoppingCriterion,
        Some(1e-5),
        11,
        code,
        None,
    )
    .unwrap();
    let key = ComputationKey::try_new(
        "nan-evidence",
        vec![sample_input(9)],
        BTreeMap::new(),
        &policy,
    )
    .unwrap();
    let artifact = b"nan-gauntlet-artifact";
    let corrupt = OutputObservation {
        artifact_hash: fs_recompute::artifact_content_hash(artifact),
        achieved_error: Some(f64::NAN),
        wall_time_s: None,
        peak_memory_bytes: None,
    };
    let mut checker = IndependentChecker::new();
    let verdict = checker.check_observation(&key, &corrupt, artifact);
    assert!(
        matches!(verdict, CheckerDisposition::InvalidEvidence { .. }),
        "NaN achieved_error must refuse as InvalidEvidence, got {verdict:?}"
    );
}
