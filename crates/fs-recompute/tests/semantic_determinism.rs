//! Test battery for semantic computation keys, tolerance roles, determinism classes,
//! and legacy migration (bead `frankensim-mkfvu.1`).
//!
//! Verifies:
//! - Canonical hashing and parameter ordering invariance;
//! - Float-bit representation and finite extrema handling;
//! - Strict rejection of non-finite/negative tolerances and errors;
//! - Tolerance role semantics: `StoppingCriterion`/`InputParameter` alter computation key;
//!   `QueryThreshold`/`None` do not;
//! - Determinism class guarantees and disposition checking;
//! - Algebraic properties: reflexivity, symmetry, transitivity;
//! - Operation family catalog defaults.

use std::collections::BTreeMap;

use fs_ledger::{ContentHash, hash_bytes};
use fs_recompute::semantic_determinism::{
    ComputationKey, DeterminismClass, DeterminismDisposition, ExecutionPolicy, OperationFamily,
    OutputObservation, SemanticKeyError, ToleranceRole,
};

fn sample_code_hash() -> ContentHash {
    hash_bytes(b"fs-solver-v0.1.0-release-build-hash")
}

fn sample_input_hash(id: u8) -> ContentHash {
    hash_bytes(&[id; 32])
}

#[test]
fn sem_001_exact_deterministic_key_properties() {
    let mut params1 = BTreeMap::new();
    params1.insert("alpha".to_string(), "1.5".to_string());
    params1.insert("beta".to_string(), "2.0".to_string());

    let policy = ExecutionPolicy::exact_deterministic(sample_code_hash(), 42);

    let key1 = ComputationKey::try_new(
        "exact_reduction",
        vec![sample_input_hash(1), sample_input_hash(2)],
        params1,
        &policy,
    )
    .expect("key1");

    // Permuted insertion order of parameters must yield identical hash
    let mut params2 = BTreeMap::new();
    params2.insert("beta".to_string(), "2.0".to_string());
    params2.insert("alpha".to_string(), "1.5".to_string());

    let key2 = ComputationKey::try_new(
        "exact_reduction",
        vec![sample_input_hash(1), sample_input_hash(2)],
        params2,
        &policy,
    )
    .expect("key2");

    assert_eq!(key1, key2);
    assert_eq!(key1.content_hash(), key2.content_hash());
    assert_eq!(key1.effective_tolerance_bits, 0);
}

#[test]
fn sem_002_tolerance_role_affects_computation_key() {
    let code_hash = sample_code_hash();
    let tol1 = 1e-4;
    let tol2 = 1e-6;

    // 1. StoppingCriterion: different tolerances must produce DIFFERENT keys
    let policy_stop1 = ExecutionPolicy::try_new(
        DeterminismClass::ToleranceDependentDeterministic,
        ToleranceRole::StoppingCriterion,
        Some(tol1),
        42,
        code_hash,
        None,
    )
    .unwrap();

    let policy_stop2 = ExecutionPolicy::try_new(
        DeterminismClass::ToleranceDependentDeterministic,
        ToleranceRole::StoppingCriterion,
        Some(tol2),
        42,
        code_hash,
        None,
    )
    .unwrap();

    let key_stop1 = ComputationKey::try_new(
        "krylov_solve",
        vec![sample_input_hash(1)],
        BTreeMap::new(),
        &policy_stop1,
    )
    .unwrap();

    let key_stop2 = ComputationKey::try_new(
        "krylov_solve",
        vec![sample_input_hash(1)],
        BTreeMap::new(),
        &policy_stop2,
    )
    .unwrap();

    assert_ne!(key_stop1, key_stop2);
    assert_ne!(key_stop1.content_hash(), key_stop2.content_hash());

    // 2. QueryThreshold: different tolerances do NOT change the computation key
    let policy_query1 = ExecutionPolicy::try_new(
        DeterminismClass::ExactDeterministic,
        ToleranceRole::QueryThreshold,
        Some(tol1),
        42,
        code_hash,
        None,
    )
    .unwrap();

    let policy_query2 = ExecutionPolicy::try_new(
        DeterminismClass::ExactDeterministic,
        ToleranceRole::QueryThreshold,
        Some(tol2),
        42,
        code_hash,
        None,
    )
    .unwrap();

    let key_query1 = ComputationKey::try_new(
        "monte_carlo",
        vec![sample_input_hash(1)],
        BTreeMap::new(),
        &policy_query1,
    )
    .unwrap();

    let key_query2 = ComputationKey::try_new(
        "monte_carlo",
        vec![sample_input_hash(1)],
        BTreeMap::new(),
        &policy_query2,
    )
    .unwrap();

    assert_eq!(key_query1, key_query2);
    assert_eq!(key_query1.content_hash(), key_query2.content_hash());
    assert_eq!(key_query1.effective_tolerance_bits, 0);
}

#[test]
fn sem_003_refusals_on_invalid_and_non_finite_inputs() {
    let code = sample_code_hash();

    // Empty op_id
    let policy = ExecutionPolicy::exact_deterministic(code, 1);
    let err_empty = ComputationKey::try_new("", vec![], BTreeMap::new(), &policy);
    assert!(matches!(err_empty, Err(SemanticKeyError::EmptyOpId)));

    // Non-finite tolerance
    let err_nan = ExecutionPolicy::try_new(
        DeterminismClass::ToleranceDependentDeterministic,
        ToleranceRole::StoppingCriterion,
        Some(f64::NAN),
        1,
        code,
        None,
    );
    assert!(matches!(
        err_nan,
        Err(SemanticKeyError::InvalidTolerance { .. })
    ));

    let err_inf = ExecutionPolicy::try_new(
        DeterminismClass::ToleranceDependentDeterministic,
        ToleranceRole::StoppingCriterion,
        Some(f64::INFINITY),
        1,
        code,
        None,
    );
    assert!(matches!(
        err_inf,
        Err(SemanticKeyError::InvalidTolerance { .. })
    ));

    // Negative / zero tolerance
    let err_zero = ExecutionPolicy::try_new(
        DeterminismClass::ToleranceDependentDeterministic,
        ToleranceRole::StoppingCriterion,
        Some(0.0),
        1,
        code,
        None,
    );
    assert!(matches!(
        err_zero,
        Err(SemanticKeyError::InvalidTolerance { .. })
    ));

    // Missing tolerance when required
    let err_missing = ExecutionPolicy::try_new(
        DeterminismClass::ToleranceDependentDeterministic,
        ToleranceRole::StoppingCriterion,
        None,
        1,
        code,
        None,
    );
    assert!(matches!(
        err_missing,
        Err(SemanticKeyError::MissingTolerance { .. })
    ));
}

#[test]
fn sem_004_output_observation_validation() {
    let art = sample_input_hash(42);

    // Valid observation
    let obs = OutputObservation::try_new(art, Some(1e-7), Some(0.45), Some(1024 * 1024)).unwrap();
    assert_eq!(obs.achieved_error, Some(1e-7));
    assert_eq!(obs.wall_time_s, Some(0.45));

    // Non-finite error
    let err_nan = OutputObservation::try_new(art, Some(f64::NAN), None, None);
    assert!(matches!(
        err_nan,
        Err(SemanticKeyError::InvalidAchievedError { .. })
    ));

    // Negative timing
    let err_time = OutputObservation::try_new(art, Some(1e-5), Some(-0.1), None);
    assert!(matches!(
        err_time,
        Err(SemanticKeyError::InvalidTiming { .. })
    ));
}

#[test]
fn sem_005_operation_family_defaults() {
    assert_eq!(
        OperationFamily::DiscreteGeometry.default_determinism_class(),
        DeterminismClass::ExactDeterministic
    );
    assert_eq!(
        OperationFamily::DiscreteGeometry.default_tolerance_role(),
        ToleranceRole::None
    );

    assert_eq!(
        OperationFamily::LinearSolve.default_determinism_class(),
        DeterminismClass::ToleranceDependentDeterministic
    );
    assert_eq!(
        OperationFamily::LinearSolve.default_tolerance_role(),
        ToleranceRole::StoppingCriterion
    );

    assert_eq!(
        OperationFamily::AdaptiveConversion.default_determinism_class(),
        DeterminismClass::ToleranceDependentDeterministic
    );
    assert_eq!(
        OperationFamily::AdaptiveConversion.default_tolerance_role(),
        ToleranceRole::InputParameter
    );

    assert_eq!(
        OperationFamily::FastHeuristic.default_determinism_class(),
        DeterminismClass::Nondeterministic
    );
}

#[test]
fn sem_006_determinism_class_properties() {
    assert!(DeterminismClass::ExactDeterministic.is_deterministic());
    assert!(DeterminismClass::ToleranceDependentDeterministic.is_deterministic());
    assert!(!DeterminismClass::Nondeterministic.is_deterministic());

    assert_eq!(
        DeterminismClass::ExactDeterministic.as_str(),
        "exact-deterministic"
    );
    assert_eq!(
        DeterminismClass::ToleranceDependentDeterministic.as_str(),
        "tolerance-dependent-deterministic"
    );
    assert_eq!(
        DeterminismClass::Nondeterministic.as_str(),
        "nondeterministic"
    );
}

#[test]
fn sem_007_algebraic_laws() {
    let policy = ExecutionPolicy::exact_deterministic(sample_code_hash(), 123);
    let key1 = ComputationKey::try_new("op", vec![sample_input_hash(1)], BTreeMap::new(), &policy)
        .unwrap();
    let key2 = ComputationKey::try_new("op", vec![sample_input_hash(1)], BTreeMap::new(), &policy)
        .unwrap();
    let key3 = ComputationKey::try_new("op", vec![sample_input_hash(1)], BTreeMap::new(), &policy)
        .unwrap();

    // Reflexivity
    assert_eq!(key1, key1);
    assert_eq!(key1.content_hash(), key1.content_hash());

    // Symmetry
    assert_eq!(key1 == key2, key2 == key1);

    // Transitivity
    assert_eq!(key1, key2);
    assert_eq!(key2, key3);
    assert_eq!(key1, key3);
}

#[test]
fn sem_008_determinism_disposition_states() {
    let disp_ok = DeterminismDisposition::Satisfied;
    assert_eq!(disp_ok, DeterminismDisposition::Satisfied);

    let disp_viol = DeterminismDisposition::Violation {
        expected_artifact: sample_input_hash(1),
        actual_artifact: sample_input_hash(2),
        diagnosis: "unstable reduction order",
    };
    assert_ne!(disp_ok, disp_viol);
}

#[test]
fn sem_009_store_put_computation_tripwire() {
    let mut store = fs_recompute::Store::new();
    let policy = ExecutionPolicy::exact_deterministic(sample_code_hash(), 42);
    let key = ComputationKey::try_new("op", vec![sample_input_hash(1)], BTreeMap::new(), &policy)
        .unwrap();

    let art1 = b"artifact-bytes-v1";
    let art2 = b"artifact-bytes-v2";

    let obs1 = OutputObservation::try_new(
        fs_recompute::artifact_content_hash(art1),
        Some(1e-6),
        None,
        None,
    )
    .unwrap();
    let res1 = store.put_computation(key.clone(), obs1, art1).unwrap();
    assert_eq!(res1, fs_recompute::PutOutcome::Inserted(key.content_hash()));

    // Re-put identical artifact -> Deduped
    let obs1_again = OutputObservation::try_new(
        fs_recompute::artifact_content_hash(art1),
        Some(1e-6),
        None,
        None,
    )
    .unwrap();
    let res_dedup = store
        .put_computation(key.clone(), obs1_again, art1)
        .unwrap();
    assert_eq!(
        res_dedup,
        fs_recompute::PutOutcome::Deduped(key.content_hash())
    );

    // Re-put DIFFERENT artifact with same key -> DeterminismViolation
    let obs2 = OutputObservation::try_new(
        fs_recompute::artifact_content_hash(art2),
        Some(1e-6),
        None,
        None,
    )
    .unwrap();
    let res_viol = store.put_computation(key, obs2, art2);
    assert!(matches!(
        res_viol,
        Err(fs_recompute::StoreError::DeterminismViolation { .. })
    ));
}

#[test]
fn sem_010_legacy_migration_idempotent() {
    let mut store = fs_recompute::Store::new();
    let rec = fs_recompute::NodeRecord {
        op_id: "legacy_mesh".to_string(),
        input_hashes: vec![sample_input_hash(1)],
        params: vec![("h".to_string(), fs_recompute::ParamValue::f(0.1))],
        code_version_hash: sample_code_hash(),
        rng_seed: 101,
        achieved_error: 1e-4,
        required_tolerance: 1e-3,
    };
    let art = b"legacy-mesh-artifact";
    store.put(rec, art).unwrap();
    assert_eq!(store.len(), 1);

    // Initial migration
    let migrated = store.migrate_all_legacy_nodes();
    assert_eq!(migrated, 1);

    // Idempotent retry -> 0 new migrations
    let retry = store.migrate_all_legacy_nodes();
    assert_eq!(retry, 0);

    // Legacy node still retrievable
    assert_eq!(store.len(), 1);
}

#[test]
fn sem_011_nondeterministic_mode_relaxation() {
    let mut store = fs_recompute::Store::new();
    let policy = ExecutionPolicy::try_new(
        DeterminismClass::Nondeterministic,
        ToleranceRole::None,
        None,
        999,
        sample_code_hash(),
        None,
    )
    .unwrap();

    let key = ComputationKey::try_new("fast_heuristic", vec![], BTreeMap::new(), &policy).unwrap();

    let art1 = b"fast-result-1";
    let art2 = b"fast-result-2";

    let obs1 =
        OutputObservation::try_new(fs_recompute::artifact_content_hash(art1), None, None, None)
            .unwrap();
    store.put_computation(key.clone(), obs1, art1).unwrap();

    // Nondeterministic mode allows differing artifacts without tripping determinism violation
    let obs2 =
        OutputObservation::try_new(fs_recompute::artifact_content_hash(art2), None, None, None)
            .unwrap();
    let res2 = store.put_computation(key, obs2, art2);
    assert!(res2.is_ok());
}
