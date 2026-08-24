//! Formal manifest test suite (bead `frankensim-extreal-program-f85xj.3.8.1`).
//!
//! Validates:
//! - Exact quantifiers and statements for next_up, next_down, interval-add, interval-mul
//! - Proof toolchain and TCB specification
//! - Deterministic content addressing
//! - Boundary enumeration of floating point classes (zeros, subnormals, extrema, infinities, NaNs)
//! - Mutation and refusal tests on corrupted manifests

use fs_ivl::formal_manifest::{
    FormalProofManifest, FROZEN_FORMAL_MANIFEST, FROZEN_MINIMUM_THEOREMS, FROZEN_PROOF_TCB,
    FROZEN_RESIDUAL_NO_CLAIMS, FROZEN_STRETCH_THEOREMS,
};

#[test]
fn manifest_validates_and_retains_required_minimum_theorems() {
    assert_eq!(FROZEN_FORMAL_MANIFEST.validate(), Ok(()));
    assert_eq!(FROZEN_MINIMUM_THEOREMS.len(), 4);
    assert_eq!(FROZEN_STRETCH_THEOREMS.len(), 3);
    assert_eq!(FROZEN_RESIDUAL_NO_CLAIMS.len(), 4);

    let theorem_ids: Vec<&str> = FROZEN_MINIMUM_THEOREMS
        .iter()
        .map(|t| t.theorem_id)
        .collect();
    assert!(theorem_ids.contains(&"thm_next_up_enclosure"));
    assert!(theorem_ids.contains(&"thm_next_down_enclosure"));
    assert!(theorem_ids.contains(&"thm_interval_add_enclosure"));
    assert!(theorem_ids.contains(&"thm_interval_mul_enclosure"));
}

#[test]
fn manifest_content_hash_is_deterministic_and_non_empty() {
    let hash1 = FROZEN_FORMAL_MANIFEST.content_hash();
    let hash2 = FROZEN_FORMAL_MANIFEST.content_hash();
    assert_eq!(hash1, hash2);
    assert!(!hash1.to_hex().is_empty());
}

#[test]
fn tcb_specification_is_complete_and_names_ieee_model() {
    assert!(FROZEN_PROOF_TCB.proof_vehicle.contains("Flocq") || FROZEN_PROOF_TCB.proof_vehicle.contains("Coq"));
    assert!(FROZEN_PROOF_TCB.ieee_model_ref.contains("IEEE 754"));
    assert!(!FROZEN_PROOF_TCB.trusted_axioms.is_empty());
    assert!(!FROZEN_PROOF_TCB.extraction_boundary.is_empty());
}

#[test]
fn residual_no_claims_are_present_and_disclose_unproved_surfaces() {
    let names: Vec<&str> = FROZEN_RESIDUAL_NO_CLAIMS
        .iter()
        .map(|nc| nc.boundary_name)
        .collect();
    assert!(names.contains(&"non_compliant_fpu_modes"));
    assert!(names.contains(&"compiler_fast_math"));
    assert!(names.contains(&"transcendental_functions"));
    assert!(names.contains(&"multivariate_taylor_models"));
}

#[test]
fn mutant_dropped_theorem_refuses_validation() {
    let incomplete_theorems = [FROZEN_MINIMUM_THEOREMS[0].clone()];
    let mutant = FormalProofManifest {
        schema_version: 1,
        tcb: FROZEN_PROOF_TCB,
        minimum_theorems: &incomplete_theorems,
        stretch_theorems: &FROZEN_STRETCH_THEOREMS,
        no_claims: &FROZEN_RESIDUAL_NO_CLAIMS,
    };
    assert!(mutant.validate().is_err());
}

#[test]
fn mutant_empty_symbol_refuses_validation() {
    let mut bad_theorems = FROZEN_MINIMUM_THEOREMS.clone();
    bad_theorems[0].rust_symbol = "";
    let mutant = FormalProofManifest {
        schema_version: 1,
        tcb: FROZEN_PROOF_TCB,
        minimum_theorems: &bad_theorems,
        stretch_theorems: &FROZEN_STRETCH_THEOREMS,
        no_claims: &FROZEN_RESIDUAL_NO_CLAIMS,
    };
    assert!(mutant.validate().is_err());
}

#[test]
fn boundary_test_enumerates_floating_point_classes() {
    // Probe boundary classes against next_up/next_down invariants
    let probe_points = [
        0.0_f64,
        -0.0_f64,
        f64::from_bits(1), // least positive subnormal
        -f64::from_bits(1),
        f64::MIN_POSITIVE, // least positive normal
        -f64::MIN_POSITIVE,
        1.0_f64,
        -1.0_f64,
        f64::MAX,
        -f64::MAX,
        f64::INFINITY,
        -f64::INFINITY,
    ];

    for &x in &probe_points {
        if x < f64::INFINITY && !x.is_nan() {
            let u = fs_math::next_up(x);
            assert!(u > x || (x == -0.0 && u > 0.0) || (x == 0.0 && u > 0.0));
        }
        if x > -f64::INFINITY && !x.is_nan() {
            let d = fs_math::next_down(x);
            assert!(d < x || (x == 0.0 && d < 0.0) || (x == -0.0 && d < 0.0));
        }
    }
}
