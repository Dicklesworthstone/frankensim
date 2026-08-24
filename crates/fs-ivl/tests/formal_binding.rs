//! Formal binding and deliberate divergence test battery
//! (bead `frankensim-extreal-program-f85xj.3.8.3`).
//!
//! Validates:
//! - Exact symbol and file bindings for the 4 core primitives
//! - Equivalence against bit-level ground truth model across f64 boundary classes
//! - Enclosure verification for interval addition and multiplication
//! - Deliberate divergence detection: injected mutants fail and produce witnesses

use fs_ivl::formal_binding::{
    binding_manifest_fingerprint, verify_boundary_class, verify_interval_add_enclosure,
    verify_interval_mul_enclosure, FROZEN_MODEL_BINDINGS,
};
use fs_ivl::interval::Interval;

#[test]
fn model_bindings_inventory_is_complete_and_deterministic() {
    assert_eq!(FROZEN_MODEL_BINDINGS.len(), 4);
    let fp1 = binding_manifest_fingerprint();
    let fp2 = binding_manifest_fingerprint();
    assert_eq!(fp1, fp2);
    assert!(!fp1.to_hex().is_empty());
}

#[test]
fn f64_boundary_classes_match_formal_model_exactly() {
    // 1. Zeros
    verify_boundary_class("zeros", &[0.0, -0.0]).expect("zeros match model");

    // 2. Subnormals
    let subnormals = [
        f64::from_bits(1),
        f64::from_bits(2),
        f64::from_bits(100),
        -f64::from_bits(1),
        -f64::from_bits(2),
    ];
    verify_boundary_class("subnormals", &subnormals).expect("subnormals match model");

    // 3. Normal-subnormal transition
    let transitions = [
        f64::MIN_POSITIVE,
        -f64::MIN_POSITIVE,
        f64::from_bits(0x000fffffffffffff), // largest subnormal
        -f64::from_bits(0x000fffffffffffff),
    ];
    verify_boundary_class("transitions", &transitions).expect("transitions match model");

    // 4. Powers of two and normals
    let normals = [1.0, -1.0, 2.0, 0.5, 1024.0, 1e-100, 1e100];
    verify_boundary_class("normals", &normals).expect("normals match model");

    // 5. Extrema and infinities
    let extrema = [f64::MAX, -f64::MAX, f64::INFINITY, -f64::INFINITY];
    verify_boundary_class("extrema", &extrema).expect("extrema match model");
}

#[test]
fn interval_add_enclosure_verified_on_test_pairs() {
    let pairs = [
        (
            Interval::new(1.0, 2.0),
            Interval::new(3.0, 4.0),
            1.5,
            3.5,
        ),
        (
            Interval::new(-1.0, 1.0),
            Interval::new(-2.0, 2.0),
            0.0,
            0.0,
        ),
        (
            Interval::new(f64::MIN_POSITIVE, 1.0),
            Interval::new(f64::MIN_POSITIVE, 2.0),
            f64::MIN_POSITIVE,
            1.0,
        ),
    ];
    assert!(verify_interval_add_enclosure(&pairs).is_ok());
}

#[test]
fn interval_mul_enclosure_verified_on_test_pairs() {
    let pairs = [
        (
            Interval::new(1.0, 2.0),
            Interval::new(3.0, 4.0),
            1.5,
            3.5,
        ),
        (
            Interval::new(-2.0, 3.0),
            Interval::new(1.0, 5.0),
            -1.0,
            2.0,
        ),
        (
            Interval::new(0.0, 1.0),
            Interval::new(0.0, 100.0),
            0.5,
            50.0,
        ),
    ];
    assert!(verify_interval_mul_enclosure(&pairs).is_ok());
}

#[test]
fn deliberate_divergence_caught_on_bad_interval_addition() {
    let outside_pair = [(
        Interval::new(1.0, 2.0),
        Interval::new(1.0, 2.0),
        3.0,
        3.0, // true sum = 6.0 is outside [2, 4]
    )];
    assert!(verify_interval_add_enclosure(&outside_pair).is_err());
}

#[test]
fn deliberate_divergence_caught_on_bad_interval_multiplication() {
    let outside_pair = [(
        Interval::new(1.0, 2.0),
        Interval::new(1.0, 2.0),
        3.0,
        3.0, // true product = 9.0 is outside [1, 4]
    )];
    assert!(verify_interval_mul_enclosure(&outside_pair).is_err());
}
