//! G0/G3/G5 coverage for canonical correlated-stack evaluation logs.

use std::num::NonZeroU64;

use fs_blake3::hash_domain;
use fs_evidence::COLOR_ALGEBRA_VERSION;
use fs_toleralloc::{
    AdmittedCorrelationModel, CORRELATED_STACK_EVALUATION_ALGORITHM_V1,
    CORRELATED_STACK_EVALUATION_IDENTITY_DOMAIN_V1, CORRELATED_STACK_EVALUATION_LOG_SCHEMA_V1,
    ColorRank, CorrelatedStackError, CorrelatedStackEvaluationLogErrorV1,
    CorrelatedStackEvaluationLogV1, CorrelatedStackTerm,
    MAX_CORRELATED_STACK_EVALUATION_LOG_BYTES_V1, MAX_CORRELATED_STACK_TERM_NAME_BYTES_V1,
    MAX_CORRELATED_STACK_TERMS_V1, ScalarIssue, propagate_correlated_stack_logged,
};

fn model(namespace: &str, version: u64, digest_byte: u8, rho: f64) -> AdmittedCorrelationModel {
    AdmittedCorrelationModel::try_new(
        namespace,
        NonZeroU64::new(version).expect("fixture versions are nonzero"),
        [digest_byte; 32],
        2,
        vec![1.0, 0.0, rho, 0.6],
    )
    .expect("manufactured two-axis model is admissible")
}

fn term(
    name: &str,
    signed_sensitivity: f64,
    sensitivity_color: ColorRank,
    standard_deviation: f64,
) -> CorrelatedStackTerm {
    CorrelatedStackTerm {
        name: name.to_string(),
        signed_sensitivity,
        sensitivity_color,
        standard_deviation,
    }
}

fn base_terms() -> Vec<CorrelatedStackTerm> {
    vec![
        term("pitch-error", 1.25, ColorRank::Verified, 0.75),
        term("radial-runout", -0.5, ColorRank::Estimated, 2.0),
    ]
}

fn base_log() -> CorrelatedStackEvaluationLogV1 {
    propagate_correlated_stack_logged(&model("gear/evaluation-log", 7, 0x5a, 0.8), &base_terms())
        .expect("base evaluation logs")
}

fn push_u32(output: &mut Vec<u8>, value: u32) {
    output.extend_from_slice(&value.to_le_bytes());
}

fn push_u64(output: &mut Vec<u8>, value: u64) {
    output.extend_from_slice(&value.to_le_bytes());
}

fn push_len_bytes(output: &mut Vec<u8>, value: &[u8]) {
    push_u64(
        output,
        u64::try_from(value.len()).expect("fixture lengths fit u64"),
    );
    output.extend_from_slice(value);
}

fn push_f64(output: &mut Vec<u8>, value: f64) {
    push_u64(output, value.to_bits());
}

fn independent_color_tag(rank: ColorRank) -> u8 {
    match rank {
        ColorRank::Estimated => 1,
        ColorRank::Validated => 2,
        ColorRank::Verified => 3,
    }
}

/// Independent specification oracle: this intentionally does not call any
/// production encoder or size helper.
fn independent_preimage(log: &CorrelatedStackEvaluationLogV1) -> Vec<u8> {
    let mut output = Vec::new();
    output.extend_from_slice(b"FSTLOGV1");
    push_u32(&mut output, 1);
    push_u32(&mut output, 1);
    push_u32(&mut output, 2);
    push_len_bytes(
        &mut output,
        b"org.frankensim.fs-toleralloc.correlated-stack-evaluation-log.v1",
    );

    let model = log.model();
    push_len_bytes(&mut output, model.namespace().as_bytes());
    push_u64(&mut output, model.schema_version().get());
    output.extend_from_slice(&model.semantic_digest());
    push_u64(
        &mut output,
        u64::try_from(model.dimension()).expect("bounded dimension fits u64"),
    );
    push_u64(
        &mut output,
        u64::try_from(model.lower_factor().len()).expect("bounded factor fits u64"),
    );
    for &factor in model.lower_factor() {
        push_f64(&mut output, factor);
    }
    push_f64(&mut output, model.max_row_norm_defect());

    push_u64(
        &mut output,
        u64::try_from(log.terms().len()).expect("bounded terms fit u64"),
    );
    for (ordinal, term) in log.terms().iter().enumerate() {
        push_u64(
            &mut output,
            u64::try_from(ordinal).expect("bounded ordinal fits u64"),
        );
        push_len_bytes(&mut output, term.name.as_bytes());
        push_f64(&mut output, term.signed_sensitivity);
        output.push(independent_color_tag(term.sensitivity_color));
        push_f64(&mut output, term.standard_deviation);
    }
    push_f64(&mut output, log.independent_standard_deviation());
    push_f64(&mut output, log.independent_variance());
    push_f64(&mut output, log.correlated_standard_deviation());
    push_f64(&mut output, log.correlated_variance());
    push_f64(&mut output, log.correlation_variance_delta());
    output
}

#[test]
fn g0_complete_nontrivial_preimage_matches_an_independent_oracle() {
    let log = base_log();
    let expected = independent_preimage(&log);
    assert_eq!(CORRELATED_STACK_EVALUATION_LOG_SCHEMA_V1, 1);
    assert_eq!(CORRELATED_STACK_EVALUATION_ALGORITHM_V1, 1);
    assert_eq!(COLOR_ALGEBRA_VERSION, 2);
    assert_eq!(
        CORRELATED_STACK_EVALUATION_IDENTITY_DOMAIN_V1,
        "org.frankensim.fs-toleralloc.correlated-stack-evaluation-log.v1",
    );
    assert_eq!(log.canonical_preimage(), expected);

    let expected_identity = hash_domain(
        "org.frankensim.fs-toleralloc.correlated-stack-evaluation-log.v1",
        &expected,
    );
    let identity = log.identity();
    assert_eq!(identity.as_bytes(), expected_identity.as_bytes());
    assert_eq!(identity.to_hex(), expected_identity.to_hex());

    assert_eq!(log.receipt().model(), log.model());
    assert_eq!(log.receipt().terms(), log.terms());
    assert_eq!(
        log.receipt().independent_standard_deviation().to_bits(),
        log.independent_standard_deviation().to_bits(),
    );
    assert_eq!(
        log.receipt().independent_variance().to_bits(),
        log.independent_variance().to_bits(),
    );
    assert_eq!(
        log.receipt().correlated_standard_deviation().to_bits(),
        log.correlated_standard_deviation().to_bits(),
    );
    assert_eq!(
        log.receipt().correlated_variance().to_bits(),
        log.correlated_variance().to_bits(),
    );
    assert_eq!(
        log.receipt().correlation_variance_delta().to_bits(),
        log.correlation_variance_delta().to_bits(),
    );
}

#[test]
fn g3_every_caller_semantic_input_moves_the_candidate_identity() {
    let baseline = base_log().identity();
    let base = base_terms();

    let cases = [
        propagate_correlated_stack_logged(&model("gear/evaluation-log-v2", 7, 0x5a, 0.8), &base),
        propagate_correlated_stack_logged(&model("gear/evaluation-log", 8, 0x5a, 0.8), &base),
        propagate_correlated_stack_logged(&model("gear/evaluation-log", 7, 0x6b, 0.8), &base),
        propagate_correlated_stack_logged(&model("gear/evaluation-log", 7, 0x5a, -0.8), &base),
        propagate_correlated_stack_logged(
            &model("gear/evaluation-log", 7, 0x5a, 0.8),
            &[
                term("pitch-error-renamed", 1.25, ColorRank::Verified, 0.75),
                base[1].clone(),
            ],
        ),
        propagate_correlated_stack_logged(
            &model("gear/evaluation-log", 7, 0x5a, 0.8),
            &[
                term("pitch-error", 1.25, ColorRank::Validated, 0.75),
                base[1].clone(),
            ],
        ),
        propagate_correlated_stack_logged(
            &model("gear/evaluation-log", 7, 0x5a, 0.8),
            &[
                term("pitch-error", 1.5, ColorRank::Verified, 0.75),
                base[1].clone(),
            ],
        ),
        propagate_correlated_stack_logged(
            &model("gear/evaluation-log", 7, 0x5a, 0.8),
            &[
                term("pitch-error", 1.25, ColorRank::Verified, 0.5),
                base[1].clone(),
            ],
        ),
    ];
    for case in cases {
        assert_ne!(
            case.expect("semantic mutation remains valid").identity(),
            baseline
        );
    }

    let three_axis_model = AdmittedCorrelationModel::try_new(
        "gear/evaluation-log",
        NonZeroU64::new(7).expect("seven is nonzero"),
        [0x5a; 32],
        3,
        vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
    )
    .expect("three-axis identity factor is admissible");
    let three_axis = propagate_correlated_stack_logged(
        &three_axis_model,
        &[
            term("pitch-error", 1.25, ColorRank::Verified, 0.75),
            term("radial-runout", -0.5, ColorRank::Estimated, 2.0),
            term("lead-error", 0.25, ColorRank::Validated, 0.5),
        ],
    )
    .expect("dimension/count mutation remains valid");
    assert_ne!(three_axis.identity(), baseline);

    let equal_values = [
        term("axis-a", 1.0, ColorRank::Validated, 1.0),
        term("axis-b", 1.0, ColorRank::Validated, 1.0),
    ];
    let first = propagate_correlated_stack_logged(
        &model("gear/evaluation-log", 7, 0x5a, 0.8),
        &equal_values,
    )
    .expect("first order logs");
    let reversed_names = [equal_values[1].clone(), equal_values[0].clone()];
    let second = propagate_correlated_stack_logged(
        &model("gear/evaluation-log", 7, 0x5a, 0.8),
        &reversed_names,
    )
    .expect("second order logs");
    assert_eq!(
        first.correlated_variance().to_bits(),
        second.correlated_variance().to_bits(),
        "the positional reorder witness deliberately preserves the numeric result",
    );
    assert_ne!(first.identity(), second.identity());
}

#[test]
fn g5_retained_model_and_terms_replay_bit_exactly() {
    let first = base_log();
    let replay = propagate_correlated_stack_logged(first.model(), first.terms())
        .expect("retained inputs replay");
    assert_eq!(first.receipt(), replay.receipt());
    assert_eq!(first.canonical_preimage(), replay.canonical_preimage());
    assert_eq!(first.identity(), replay.identity());
}

fn identity_factor(dimension: usize) -> Vec<f64> {
    let mut factor = vec![0.0; dimension * dimension];
    for index in 0..dimension {
        factor[index * dimension + index] = 1.0;
    }
    factor
}

fn max_width_terms(count: usize) -> Vec<CorrelatedStackTerm> {
    (0..count)
        .map(|index| {
            let prefix = format!("{index:03}-");
            let name = format!(
                "{prefix}{}",
                "x".repeat(MAX_CORRELATED_STACK_TERM_NAME_BYTES_V1 - prefix.len()),
            );
            term(&name, 1.0, ColorRank::Verified, 1.0)
        })
        .collect()
}

#[test]
fn g0_exact_128_axis_envelope_and_invalid_inputs_publish_no_log() {
    let namespace = format!(
        "gear/{}",
        "a".repeat(fs_toleralloc::MAX_CORRELATION_MODEL_NAMESPACE_BYTES_V1 - "gear/".len()),
    );
    let max_model = AdmittedCorrelationModel::try_new(
        namespace,
        NonZeroU64::new(1).expect("one is nonzero"),
        [0x7f; 32],
        MAX_CORRELATED_STACK_TERMS_V1,
        identity_factor(MAX_CORRELATED_STACK_TERMS_V1),
    )
    .expect("maximum identity factor is admissible");
    let max_terms = max_width_terms(MAX_CORRELATED_STACK_TERMS_V1);
    let max_log = propagate_correlated_stack_logged(&max_model, &max_terms)
        .expect("exact maximum envelope logs");
    assert_eq!(max_log.terms().len(), MAX_CORRELATED_STACK_TERMS_V1);
    assert_eq!(
        max_log.canonical_preimage().len(),
        MAX_CORRELATED_STACK_EVALUATION_LOG_BYTES_V1,
    );

    let too_many = max_width_terms(MAX_CORRELATED_STACK_TERMS_V1 + 1);
    assert_eq!(
        propagate_correlated_stack_logged(&max_model, &too_many),
        Err(CorrelatedStackEvaluationLogErrorV1::Stack(
            CorrelatedStackError::TooManyTerms {
                actual: MAX_CORRELATED_STACK_TERMS_V1 + 1,
                max: MAX_CORRELATED_STACK_TERMS_V1,
            },
        )),
    );
    assert_eq!(
        propagate_correlated_stack_logged(&max_model, &[]),
        Err(CorrelatedStackEvaluationLogErrorV1::Stack(
            CorrelatedStackError::NoTerms,
        )),
    );
    assert_eq!(
        propagate_correlated_stack_logged(&max_model, &max_terms[..1]),
        Err(CorrelatedStackEvaluationLogErrorV1::Stack(
            CorrelatedStackError::DimensionMismatch {
                model: MAX_CORRELATED_STACK_TERMS_V1,
                terms: 1,
            },
        )),
    );

    let two_axis = model("gear/evaluation-log", 1, 0x31, 0.8);
    let duplicate = [
        term("same", 1.0, ColorRank::Verified, 1.0),
        term("SAME", 1.0, ColorRank::Verified, 1.0),
    ];
    assert!(matches!(
        propagate_correlated_stack_logged(&two_axis, &duplicate),
        Err(CorrelatedStackEvaluationLogErrorV1::Stack(
            CorrelatedStackError::AmbiguousTermName { .. }
        )),
    ));
    let nonfinite = [
        term("a", f64::NAN, ColorRank::Verified, 1.0),
        term("b", 1.0, ColorRank::Verified, 1.0),
    ];
    assert!(matches!(
        propagate_correlated_stack_logged(&two_axis, &nonfinite),
        Err(CorrelatedStackEvaluationLogErrorV1::Stack(
            CorrelatedStackError::InvalidTermField {
                index: 0,
                field: "signed_sensitivity",
                issue: ScalarIssue::NonFinite,
                ..
            }
        )),
    ));
}

#[test]
fn g0_exact_positive_zero_stack_is_logged_and_negative_zero_is_refused() {
    let model = model("gear/evaluation-log-zero", 1, 0x42, 0.8);
    let positive_zero = [
        term("a", 0.0, ColorRank::Estimated, 1.0),
        term("b", 0.0, ColorRank::Verified, 2.0),
    ];
    let log = propagate_correlated_stack_logged(&model, &positive_zero)
        .expect("all-positive-zero sensitivities remain a bound evaluation");
    for value in [
        log.independent_standard_deviation(),
        log.independent_variance(),
        log.correlated_standard_deviation(),
        log.correlated_variance(),
        log.correlation_variance_delta(),
    ] {
        assert_eq!(value.to_bits(), 0.0_f64.to_bits());
    }
    assert_eq!(log.model(), &model);
    assert_eq!(log.terms(), positive_zero);
    assert_eq!(log.canonical_preimage(), independent_preimage(&log));
    let replay = propagate_correlated_stack_logged(log.model(), log.terms())
        .expect("retained positive-zero evaluation replays");
    assert_eq!(log.receipt(), replay.receipt());
    assert_eq!(log.canonical_preimage(), replay.canonical_preimage());
    assert_eq!(log.identity(), replay.identity());

    let negative_zero = [
        term("a", -0.0, ColorRank::Estimated, 1.0),
        term("b", 0.0, ColorRank::Verified, 2.0),
    ];
    assert!(matches!(
        propagate_correlated_stack_logged(&model, &negative_zero),
        Err(CorrelatedStackEvaluationLogErrorV1::Stack(
            CorrelatedStackError::InvalidTermField {
                index: 0,
                field: "signed_sensitivity",
                issue: ScalarIssue::NonCanonicalNegativeZero,
                ..
            }
        )),
    ));
}

#[test]
fn g3_domain_separation_is_not_plain_preimage_hashing() {
    let log = base_log();
    let domain_identity = hash_domain(
        CORRELATED_STACK_EVALUATION_IDENTITY_DOMAIN_V1,
        log.canonical_preimage(),
    );
    let neighboring_identity = hash_domain(
        "org.frankensim.fs-toleralloc.correlated-stack-evaluation-log-neighbor.v1",
        log.canonical_preimage(),
    );
    assert_eq!(log.identity().as_bytes(), domain_identity.as_bytes());
    assert_ne!(log.identity().as_bytes(), neighboring_identity.as_bytes());
}

#[test]
fn g3_each_versioned_header_component_moves_the_candidate_digest() {
    let log = base_log();
    let baseline = log.identity();
    let mut mutation_offsets = vec![0, 8, 12, 16, 20, 28];
    mutation_offsets.push(log.canonical_preimage().len() - 1);
    for offset in mutation_offsets {
        let mut moved = log.canonical_preimage().to_vec();
        moved[offset] ^= 1;
        let moved_identity = hash_domain(CORRELATED_STACK_EVALUATION_IDENTITY_DOMAIN_V1, &moved);
        assert_ne!(moved_identity.as_bytes(), baseline.as_bytes());
    }
}
