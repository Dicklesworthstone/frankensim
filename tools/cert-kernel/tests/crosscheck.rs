//! G0/G3 independent-kernel comparison and seeded discrepancy battery.

use frankensim_cert_kernel::crosscheck::run_audit;

#[test]
fn exact_and_boundary_corpora_agree_without_laundering_divergence() {
    let report = run_audit(4096).expect("deterministic corpus must be admissible");
    assert!(report.seeded_tripwire_detected);
    assert!(report.is_green(), "{}", report.json_lines());
    assert_eq!(report.operations.len(), 7);
    for operation in &report.operations {
        assert_eq!(operation.compatibility_cases, 4096);
        assert_eq!(operation.non_overlaps, 0);
        assert!(operation.exact_reference_cases > 0);
        assert_eq!(operation.exact_reference_misses, 0);
    }
    assert_eq!(
        report
            .json_lines()
            .matches("\"first_non_overlap\":null")
            .count(),
        7
    );
    println!("{}", report.json_lines());
}

#[test]
fn deterministic_artifact_is_byte_stable() {
    let first = run_audit(128).expect("first audit").json_lines();
    let second = run_audit(128).expect("second audit").json_lines();
    assert_eq!(first, second);
}

#[test]
fn zero_sample_request_refuses_before_emitting_a_report() {
    let error = run_audit(0).expect_err("zero work must refuse");
    assert_eq!(error.to_string(), "samples_per_operation must be positive");
}
