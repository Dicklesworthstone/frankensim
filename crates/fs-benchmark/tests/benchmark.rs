//! G0/G3 battery for the wedge benchmark corpus. The tests pin exact evidence
//! resolution, typed denominator-backed metrics, validation failure domains,
//! and mutation sensitivity of the length-framed corpus identity.

use fs_benchmark::{
    AdmissionDecision, AdmissionReceipt, AdmissionVerifier, BENCHMARK_VERSION, BenchmarkCorpus,
    COLOR_ALGEBRA_VERSION, CORPUS_IDENTITY_SCHEMA_VERSION, Color, ColorRank, ContentHash,
    DatasetKind, EvidenceDatum, EvidenceError, EvidenceRecord, EvidenceRef, EvidenceRole,
    InstrumentedProposal, MetricDefinition, MetricError, MetricEvidenceRole, MetricKind,
    MetricRequest, MetricSchema, ProposalEvaluation, ProposalEvaluator, ProposalEvidenceRef,
    ProposalRefusal, ReferenceAuthority, accept_rate, audit, audit_corpus,
    audit_corpus_with_verifier, benchmark_corpus, conflict_rate, corpus_digest, corpus_digest_for,
    design_tasks, edit_skip_metric, edit_traces, evaluate_proposal, instrumented_proposals,
    merge_conflict_metric, merge_trials, mms_battery, query_set, rate, reconstruct_metric,
    resolve_evidence, resolve_query_reference, resolve_query_reference_with_verifier,
    retained_evidence, speedup, win_rate,
};

fn assert_digest_changed(changed: &BenchmarkCorpus<'_>) {
    assert_ne!(corpus_digest_for(changed), corpus_digest());
}

fn changed_digest(mut digest: ContentHash) -> ContentHash {
    digest.0[0] ^= 1;
    digest
}

fn role(
    record: &EvidenceRecord,
    subject_id: &'static str,
    quantity: &'static str,
    units: &'static str,
) -> MetricEvidenceRole {
    MetricEvidenceRole::exact(
        EvidenceRef::exact(record.id, record.semantic_digest()),
        subject_id,
        quantity,
        units,
    )
}

fn proposal_manifest(
    proposal: &InstrumentedProposal,
    records: &[EvidenceRecord],
) -> Vec<ProposalEvidenceRef> {
    records
        .iter()
        .map(|record| ProposalEvidenceRef {
            proposal: proposal.proposal,
            subject_id: record.subject_id,
            role: record.role,
            reference: EvidenceRef::exact(record.id, record.semantic_digest()),
        })
        .collect()
}

fn assert_proposal_values(
    corpus: &BenchmarkCorpus<'_>,
    proposal: &InstrumentedProposal,
    expected: &[f64],
) {
    let ProposalEvaluation::Available { metrics } = evaluate_proposal(corpus, proposal) else {
        panic!("complete typed proposal evidence must be available");
    };
    assert_eq!(metrics.len(), expected.len());
    for (metric, expected) in metrics.iter().zip(expected) {
        assert!((metric.value - *expected).abs() < 1e-12);
    }
}

#[test]
fn every_dataset_is_populated_and_references_remain_non_authoritative_by_default() {
    assert_eq!(query_set().len(), 3);
    assert_eq!(design_tasks().len(), 3);
    assert_eq!(edit_traces().len(), 2);
    assert_eq!(mms_battery().len(), 2);
    assert_eq!(merge_trials().len(), 2);
    assert_eq!(retained_evidence().len(), 11);

    let corpus = benchmark_corpus();
    for query in query_set() {
        assert!(!query.qoi.is_empty());
        assert!(!query.units.is_empty());
        let resolved = resolve_query_reference(&corpus, query).expect("retained query evidence");
        assert_eq!(
            resolved.answer().to_bits(),
            query.reference_answer.to_bits()
        );
        assert_ne!(resolved.evidence_digest(), ContentHash([0; 32]));
        assert_eq!(resolved.authority().rank(), ColorRank::Estimated);
        assert!(!resolved.authority().is_admitted());
        assert!(matches!(
            resolved.authority(),
            ReferenceAuthority::EstimatedDeclaration(_)
        ));
    }
}

#[test]
fn raw_metric_helpers_refuse_missing_or_invalid_denominators() {
    assert!((speedup(1000.0, 400.0).expect("valid speedup") - 2.5).abs() < 1e-12);
    assert_eq!(speedup(1000.0, 0.0), Err(MetricError::ZeroDenominator));
    assert!(matches!(
        speedup(f64::NAN, 1.0),
        Err(MetricError::NonFiniteInput { role: "baseline" })
    ));
    assert!(matches!(
        speedup(1.0, -1.0),
        Err(MetricError::NegativeInput { role: "candidate" })
    ));
    assert_eq!(
        speedup(f64::MAX, f64::from_bits(1)),
        Err(MetricError::NonFiniteResult)
    );

    assert!(
        (win_rate(&[true, true, false, true, true]).expect("non-empty rate") - 0.8).abs() < 1e-12
    );
    assert_eq!(win_rate(&[]), Err(MetricError::ZeroDenominator));
    assert!((accept_rate(35, 100).expect("nonzero attempts") - 0.35).abs() < 1e-12);
    assert_eq!(rate(1, 0), Err(MetricError::ZeroDenominator));
    assert_eq!(rate(2, 1), Err(MetricError::NumeratorExceedsDenominator));
    if usize::BITS >= 64 {
        assert_eq!(
            rate((1_u64 << 53) as usize, ((1_u64 << 53) + 1) as usize),
            Err(MetricError::InexactCount {
                role: "denominator"
            })
        );
    }
}

#[test]
fn built_in_rates_reconstruct_from_typed_numerator_and_denominator_evidence() {
    let corpus = benchmark_corpus();
    let trace = &edit_traces()[0];
    let skip = reconstruct_metric(&corpus, edit_skip_metric(trace)).expect("retained skip metric");
    assert_eq!(skip.kind, MetricKind::Rate);
    assert!((skip.value - 0.8).abs() < 1e-12);
    assert_eq!(skip.numerator_evidence, trace.correct_skips_evidence);
    assert_eq!(skip.denominator_evidence, trace.total_ops_evidence);
    assert_ne!(skip.identity, ContentHash([0; 32]));
    assert_eq!(
        skip.identity,
        reconstruct_metric(&corpus, edit_skip_metric(trace))
            .expect("deterministic replay")
            .identity
    );

    let trial = &merge_trials()[0];
    let conflict = reconstruct_metric(&corpus, merge_conflict_metric(trial))
        .expect("retained conflict metric");
    assert!((conflict.value - 0.15).abs() < 1e-12);
    assert!((conflict_rate(trial).expect("nonzero merge count") - 0.15).abs() < 1e-12);
}

#[test]
fn typed_speedup_resolves_evidence_and_refuses_zero_or_mismatched_denominators() {
    const EVALUATOR: &str = "equal-accuracy-cost-ratio-v1";
    const PROPOSAL: &str = "addendum-proposal-8:planner-kill-v1";
    let baseline = EvidenceRecord::scalar(
        "metric:baseline",
        "planner-q1",
        "baseline-cost",
        "cost-unit",
        EVALUATOR,
        PROPOSAL,
        EvidenceRole::BaselineCost,
        1000.0,
    );
    let candidate = EvidenceRecord::scalar(
        "metric:candidate",
        "planner-q1",
        "candidate-cost",
        "cost-unit",
        EVALUATOR,
        PROPOSAL,
        EvidenceRole::CandidateCost,
        400.0,
    );
    let request = MetricRequest {
        schema: MetricSchema::GenericSpeedup,
        evaluator_semantics: EVALUATOR,
        proposal_semantics: PROPOSAL,
        definition: MetricDefinition::Speedup {
            baseline: role(&baseline, "planner-q1", "baseline-cost", "cost-unit"),
            candidate: role(&candidate, "planner-q1", "candidate-cost", "cost-unit"),
        },
    };
    let records = [baseline.clone(), candidate.clone()];
    let corpus = BenchmarkCorpus {
        retained_evidence: &records,
        ..benchmark_corpus()
    };
    let metric = reconstruct_metric(&corpus, request).expect("typed speedup evidence");
    assert_eq!(metric.kind, MetricKind::Speedup);
    assert!((metric.value - 2.5).abs() < 1e-12);
    assert_eq!(
        reconstruct_metric(
            &corpus,
            MetricRequest {
                schema: MetricSchema::GenericRate,
                ..request
            }
        ),
        Err(MetricError::SchemaKindMismatch)
    );

    let missing_records = [baseline.clone()];
    let missing = BenchmarkCorpus {
        retained_evidence: &missing_records,
        ..benchmark_corpus()
    };
    assert!(matches!(
        reconstruct_metric(&missing, request),
        Err(MetricError::Evidence(EvidenceError::Missing {
            id: "metric:candidate"
        }))
    ));

    let mut tampered_candidate = candidate.clone();
    tampered_candidate.datum = EvidenceDatum::Scalar(399.0);
    let tampered_records = [baseline.clone(), tampered_candidate];
    let tampered = BenchmarkCorpus {
        retained_evidence: &tampered_records,
        ..benchmark_corpus()
    };
    assert!(matches!(
        reconstruct_metric(&tampered, request),
        Err(MetricError::Evidence(EvidenceError::Tampered {
            id: "metric:candidate"
        }))
    ));

    let zero = EvidenceRecord::scalar(
        "metric:candidate",
        "planner-q1",
        "candidate-cost",
        "cost-unit",
        EVALUATOR,
        PROPOSAL,
        EvidenceRole::CandidateCost,
        0.0,
    );
    let zero_records = [baseline.clone(), zero.clone()];
    let zero_corpus = BenchmarkCorpus {
        retained_evidence: &zero_records,
        ..benchmark_corpus()
    };
    let zero_request = MetricRequest {
        definition: MetricDefinition::Speedup {
            baseline: role(&baseline, "planner-q1", "baseline-cost", "cost-unit"),
            candidate: role(&zero, "planner-q1", "candidate-cost", "cost-unit"),
        },
        ..request
    };
    assert_eq!(
        reconstruct_metric(&zero_corpus, zero_request),
        Err(MetricError::ZeroDenominator)
    );

    let wrong_units = EvidenceRecord::scalar(
        "metric:candidate",
        "planner-q1",
        "candidate-cost",
        "seconds",
        EVALUATOR,
        PROPOSAL,
        EvidenceRole::CandidateCost,
        400.0,
    );
    let wrong_unit_records = [baseline.clone(), wrong_units.clone()];
    let wrong_unit_corpus = BenchmarkCorpus {
        retained_evidence: &wrong_unit_records,
        ..benchmark_corpus()
    };
    let wrong_unit_request = MetricRequest {
        definition: MetricDefinition::Speedup {
            baseline: role(&baseline, "planner-q1", "baseline-cost", "cost-unit"),
            candidate: role(&wrong_units, "planner-q1", "candidate-cost", "cost-unit"),
        },
        ..request
    };
    assert!(matches!(
        reconstruct_metric(&wrong_unit_corpus, wrong_unit_request),
        Err(MetricError::Evidence(EvidenceError::ContextMismatch {
            field: "units",
            ..
        }))
    ));

    let wrong_baseline_units = EvidenceRecord::scalar(
        "metric:baseline",
        "planner-q1",
        "baseline-cost",
        "seconds",
        EVALUATOR,
        PROPOSAL,
        EvidenceRole::BaselineCost,
        1000.0,
    );
    let equally_wrong_records = [wrong_baseline_units.clone(), wrong_units.clone()];
    let equally_wrong_corpus = BenchmarkCorpus {
        retained_evidence: &equally_wrong_records,
        ..benchmark_corpus()
    };
    let equally_wrong_request = MetricRequest {
        definition: MetricDefinition::Speedup {
            baseline: role(
                &wrong_baseline_units,
                "planner-q1",
                "baseline-cost",
                "cost-unit",
            ),
            candidate: role(&wrong_units, "planner-q1", "candidate-cost", "cost-unit"),
        },
        ..request
    };
    assert!(matches!(
        reconstruct_metric(&equally_wrong_corpus, equally_wrong_request),
        Err(MetricError::Evidence(EvidenceError::ContextMismatch {
            field: "units",
            ..
        }))
    ));

    let swapped_request = MetricRequest {
        definition: MetricDefinition::Speedup {
            baseline: role(&candidate, "planner-q1", "candidate-cost", "cost-unit"),
            candidate: role(&baseline, "planner-q1", "baseline-cost", "cost-unit"),
        },
        ..request
    };
    assert!(matches!(
        reconstruct_metric(&corpus, swapped_request),
        Err(MetricError::RoleMismatch {
            expected: EvidenceRole::BaselineCost,
            actual: EvidenceRole::CandidateCost,
        })
    ));

    let same_reference_request = MetricRequest {
        definition: MetricDefinition::Speedup {
            baseline: role(&baseline, "planner-q1", "baseline-cost", "cost-unit"),
            candidate: role(&baseline, "planner-q1", "candidate-cost", "cost-unit"),
        },
        ..request
    };
    assert_eq!(
        reconstruct_metric(&corpus, same_reference_request),
        Err(MetricError::SameEvidence)
    );
}

#[test]
fn query_resolution_refuses_missing_duplicate_tampered_and_wrong_context_evidence() {
    let query = &query_set()[0];

    let missing_records: Vec<_> = retained_evidence()
        .iter()
        .cloned()
        .filter(|record| record.id != query.reference_evidence.id)
        .collect();
    let missing = BenchmarkCorpus {
        retained_evidence: &missing_records,
        ..benchmark_corpus()
    };
    assert!(matches!(
        resolve_query_reference(&missing, query),
        Err(EvidenceError::Missing { .. })
    ));

    let mut tampered_records = retained_evidence().to_vec();
    let EvidenceDatum::QueryReference { answer, .. } = &mut tampered_records[0].datum else {
        panic!("built-in query evidence must be a query reference");
    };
    *answer += 1.0;
    let tampered = BenchmarkCorpus {
        retained_evidence: &tampered_records,
        ..benchmark_corpus()
    };
    assert!(matches!(
        resolve_query_reference(&tampered, query),
        Err(EvidenceError::Tampered { .. })
    ));

    let resigned_replacement = EvidenceRecord::query_reference(
        query.reference_evidence.id,
        query.id,
        query.qoi,
        query.units,
        query.reference_evaluator_semantics,
        query.reference_proposal_semantics,
        query.reference_answer,
        query.tolerance,
        query.reference_cost,
        query.reference_cost_units,
        Color::Estimated {
            estimator: "resigned-replacement-v1".to_string(),
            dispersion: f64::INFINITY,
        },
        None,
    );
    let resigned_records = [resigned_replacement];
    let resigned = BenchmarkCorpus {
        retained_evidence: &resigned_records,
        ..benchmark_corpus()
    };
    assert!(matches!(
        resolve_query_reference(&resigned, query),
        Err(EvidenceError::Tampered { .. })
    ));

    let mut duplicate_records = retained_evidence().to_vec();
    duplicate_records.push(retained_evidence()[0].clone());
    let duplicate = BenchmarkCorpus {
        retained_evidence: &duplicate_records,
        ..benchmark_corpus()
    };
    assert!(matches!(
        resolve_query_reference(&duplicate, query),
        Err(EvidenceError::Duplicate { .. })
    ));

    let wrong_context = EvidenceRecord::query_reference(
        query.reference_evidence.id,
        query.id,
        "different QoI",
        query.units,
        query.reference_evaluator_semantics,
        query.reference_proposal_semantics,
        query.reference_answer,
        query.tolerance,
        query.reference_cost,
        query.reference_cost_units,
        Color::Estimated {
            estimator: "wrong-context-v1".to_string(),
            dispersion: f64::INFINITY,
        },
        None,
    );
    let records = [wrong_context.clone()];
    let context_corpus = BenchmarkCorpus {
        retained_evidence: &records,
        ..benchmark_corpus()
    };
    let mut context_query = *query;
    context_query.reference_evidence =
        EvidenceRef::exact(wrong_context.id, wrong_context.semantic_digest());
    assert!(matches!(
        resolve_query_reference(&context_corpus, &context_query),
        Err(EvidenceError::ContextMismatch {
            field: "quantity",
            ..
        })
    ));

    let corpus = benchmark_corpus();
    let mut changed = *query;
    changed.tolerance = query.tolerance * 2.0;
    assert!(matches!(
        resolve_query_reference(&corpus, &changed),
        Err(EvidenceError::ContextMismatch {
            field: "tolerance",
            ..
        })
    ));
    changed = *query;
    changed.reference_cost = query.reference_cost + 1.0;
    assert!(matches!(
        resolve_query_reference(&corpus, &changed),
        Err(EvidenceError::ContextMismatch {
            field: "reference-cost",
            ..
        })
    ));
    changed = *query;
    changed.reference_cost_units = "seconds";
    assert!(matches!(
        resolve_query_reference(&corpus, &changed),
        Err(EvidenceError::ContextMismatch {
            field: "reference-cost-units",
            ..
        })
    ));

    for (field, tampered_datum) in [
        (
            "reference-answer",
            EvidenceDatum::QueryReference {
                answer: query.reference_answer + 1.0,
                tolerance: query.tolerance,
                reference_cost: query.reference_cost,
                reference_cost_units: query.reference_cost_units,
                color: Color::Estimated {
                    estimator: "resigned-context-v1".to_string(),
                    dispersion: f64::INFINITY,
                },
                admission_receipt: None,
            },
        ),
        (
            "tolerance",
            EvidenceDatum::QueryReference {
                answer: query.reference_answer,
                tolerance: query.tolerance * 2.0,
                reference_cost: query.reference_cost,
                reference_cost_units: query.reference_cost_units,
                color: Color::Estimated {
                    estimator: "resigned-context-v1".to_string(),
                    dispersion: f64::INFINITY,
                },
                admission_receipt: None,
            },
        ),
        (
            "reference-cost",
            EvidenceDatum::QueryReference {
                answer: query.reference_answer,
                tolerance: query.tolerance,
                reference_cost: query.reference_cost + 1.0,
                reference_cost_units: query.reference_cost_units,
                color: Color::Estimated {
                    estimator: "resigned-context-v1".to_string(),
                    dispersion: f64::INFINITY,
                },
                admission_receipt: None,
            },
        ),
    ] {
        let mut record = retained_evidence()[0].clone();
        record.datum = tampered_datum;
        let mut resigned_query = *query;
        resigned_query.reference_evidence = EvidenceRef::exact(record.id, record.semantic_digest());
        let records = [record];
        let resigned_corpus = BenchmarkCorpus {
            retained_evidence: &records,
            ..benchmark_corpus()
        };
        assert!(matches!(
            resolve_query_reference(&resigned_corpus, &resigned_query),
            Err(EvidenceError::ContextMismatch {
                field: actual,
                ..
            }) if actual == field
        ));
    }
}

struct ReceiptPolicyVerifier;

impl AdmissionVerifier for ReceiptPolicyVerifier {
    fn verify(&self, _candidate: &Color, receipt: &AdmissionReceipt) -> AdmissionDecision {
        AdmissionDecision::accept(receipt.policy_fingerprint())
    }
}

#[test]
fn positive_query_color_requires_injected_admission_even_after_rehashing() {
    let color = Color::Verified {
        lo: 299.75,
        hi: 300.25,
    };
    let provisional = EvidenceRecord::query_reference(
        "query:positive",
        "positive-q1",
        "temperature",
        "K",
        "reference-evaluator-v1",
        "reference-proposal-v1",
        300.0,
        0.25,
        42.0,
        "work-unit",
        color.clone(),
        Some(AdmissionReceipt::from_parts(
            ContentHash([0; 32]),
            7,
            COLOR_ALGEBRA_VERSION,
            ContentHash([2; 32]),
        )),
    );
    let receipt = AdmissionReceipt::from_parts(
        provisional
            .query_admission_node_hash()
            .expect("positive query context identity"),
        7,
        COLOR_ALGEBRA_VERSION,
        ContentHash([2; 32]),
    );
    let record = EvidenceRecord::query_reference(
        "query:positive",
        "positive-q1",
        "temperature",
        "K",
        "reference-evaluator-v1",
        "reference-proposal-v1",
        300.0,
        0.25,
        42.0,
        "work-unit",
        color,
        Some(receipt),
    );
    let query = fs_benchmark::QueryCase {
        id: "positive-q1",
        qoi: "temperature",
        units: "K",
        tolerance: 0.25,
        reference_answer: 300.0,
        reference_cost: 42.0,
        reference_cost_units: "work-unit",
        reference_evidence: EvidenceRef::exact(record.id, record.semantic_digest()),
        reference_evaluator_semantics: "reference-evaluator-v1",
        reference_proposal_semantics: "reference-proposal-v1",
    };
    let records = [record];
    let queries = [query];
    let corpus = BenchmarkCorpus {
        query_set: &queries,
        retained_evidence: &records,
        ..benchmark_corpus()
    };

    assert!(matches!(
        resolve_query_reference(&corpus, &query),
        Err(EvidenceError::AdmissionRejected { .. })
    ));
    let admitted = resolve_query_reference_with_verifier(&corpus, &query, &ReceiptPolicyVerifier)
        .expect("injected receipt-policy verifier admits the candidate");
    assert!(admitted.authority().is_admitted());
    assert_eq!(admitted.authority().rank(), ColorRank::Verified);
    let ReferenceAuthority::Admitted(admitted_reference) = admitted.authority() else {
        unreachable!("positive query must return admitted authority");
    };
    assert_eq!(
        admitted_reference.query_context_hash(),
        admitted.query_context_hash()
    );

    let deny_all_audit = audit_corpus(&corpus);
    assert!(deny_all_audit.gaps.iter().any(|gap| {
        gap.starts_with("query 'positive-q1':") && gap.contains("was not admitted")
    }));
    let admitted_audit = audit_corpus_with_verifier(&corpus, &ReceiptPolicyVerifier);
    assert!(
        !admitted_audit
            .gaps
            .iter()
            .any(|gap| gap.starts_with("query 'positive-q1':")),
        "injected authority should admit the valid positive query: {:?}",
        admitted_audit.gaps
    );

    let mut replayed_record = records[0].clone();
    let EvidenceDatum::QueryReference { answer, .. } = &mut replayed_record.datum else {
        unreachable!("positive query record shape");
    };
    *answer = 300.1;
    let mut replayed_query = query;
    replayed_query.reference_answer = 300.1;
    replayed_query.reference_evidence =
        EvidenceRef::exact(replayed_record.id, replayed_record.semantic_digest());
    let replayed_records = [replayed_record];
    let replayed_queries = [replayed_query];
    let replayed_corpus = BenchmarkCorpus {
        query_set: &replayed_queries,
        retained_evidence: &replayed_records,
        ..benchmark_corpus()
    };
    assert!(matches!(
        resolve_query_reference_with_verifier(
            &replayed_corpus,
            &replayed_query,
            &ReceiptPolicyVerifier,
        ),
        Err(EvidenceError::AdmissionContextMismatch { .. })
    ));
}

#[test]
fn query_role_and_verified_answer_are_fail_closed() {
    let query = query_set()[0];
    let mut wrong_role = retained_evidence()[0].clone();
    wrong_role.role = EvidenceRole::BaselineCost;
    let mut resigned_query = query;
    resigned_query.reference_evidence =
        EvidenceRef::exact(wrong_role.id, wrong_role.semantic_digest());
    let wrong_role_records = [wrong_role];
    let wrong_role_queries = [resigned_query];
    let wrong_role_corpus = BenchmarkCorpus {
        query_set: &wrong_role_queries,
        retained_evidence: &wrong_role_records,
        ..benchmark_corpus()
    };
    assert!(matches!(
        resolve_query_reference(&wrong_role_corpus, &resigned_query),
        Err(EvidenceError::ContextMismatch {
            field: "role-datum",
            ..
        })
    ));

    let outside_interval = EvidenceRecord::query_reference(
        "query:outside-verified-interval",
        "outside-q1",
        "temperature",
        "K",
        "reference-evaluator-v1",
        "reference-proposal-v1",
        301.0,
        0.25,
        42.0,
        "work-unit",
        Color::Verified {
            lo: 299.75,
            hi: 300.25,
        },
        Some(AdmissionReceipt::from_parts(
            ContentHash([0; 32]),
            7,
            COLOR_ALGEBRA_VERSION,
            ContentHash([2; 32]),
        )),
    );
    let outside_interval_records = [outside_interval.clone()];
    let outside_interval_corpus = BenchmarkCorpus {
        retained_evidence: &outside_interval_records,
        ..benchmark_corpus()
    };
    assert!(matches!(
        resolve_evidence(
            &outside_interval_corpus,
            EvidenceRef::exact(outside_interval.id, outside_interval.semantic_digest(),),
        ),
        Err(EvidenceError::ContextMismatch {
            field: "verified-answer",
            ..
        })
    ));
}

#[test]
fn query_evidence_identity_binds_color_and_every_receipt_field() {
    let make = |color: Color, receipt: AdmissionReceipt| {
        EvidenceRecord::query_reference(
            "query:identity",
            "identity-q1",
            "temperature",
            "K",
            "identity-evaluator-v1",
            "identity-proposal-v1",
            300.0,
            0.25,
            42.0,
            "work-unit",
            color,
            Some(receipt),
        )
    };
    let receipt = AdmissionReceipt::from_parts(
        ContentHash([1; 32]),
        7,
        COLOR_ALGEBRA_VERSION,
        ContentHash([2; 32]),
    );
    let base = make(
        Color::Verified {
            lo: 299.75,
            hi: 300.25,
        },
        receipt,
    )
    .semantic_digest();
    assert_eq!(
        base,
        make(
            Color::Verified {
                lo: 299.75,
                hi: 300.25,
            },
            receipt,
        )
        .semantic_digest()
    );
    for changed in [
        make(
            Color::Verified {
                lo: 299.5,
                hi: 300.25,
            },
            receipt,
        ),
        make(
            Color::Verified {
                lo: 299.75,
                hi: 300.5,
            },
            receipt,
        ),
        make(
            Color::Verified {
                lo: 299.75,
                hi: 300.25,
            },
            AdmissionReceipt::from_parts(
                ContentHash([3; 32]),
                7,
                COLOR_ALGEBRA_VERSION,
                ContentHash([2; 32]),
            ),
        ),
        make(
            Color::Verified {
                lo: 299.75,
                hi: 300.25,
            },
            AdmissionReceipt::from_parts(
                ContentHash([1; 32]),
                8,
                COLOR_ALGEBRA_VERSION,
                ContentHash([2; 32]),
            ),
        ),
        make(
            Color::Verified {
                lo: 299.75,
                hi: 300.25,
            },
            AdmissionReceipt::from_parts(
                ContentHash([1; 32]),
                7,
                COLOR_ALGEBRA_VERSION + 1,
                ContentHash([2; 32]),
            ),
        ),
        make(
            Color::Verified {
                lo: 299.75,
                hi: 300.25,
            },
            AdmissionReceipt::from_parts(
                ContentHash([1; 32]),
                7,
                COLOR_ALGEBRA_VERSION,
                ContentHash([4; 32]),
            ),
        ),
    ] {
        assert_ne!(changed.semantic_digest(), base);
    }
}

#[test]
fn count_payload_tampering_and_duplicate_authoritative_digests_are_refused() {
    let reference = edit_traces()[0].correct_skips_evidence;
    let mut tampered_records = retained_evidence().to_vec();
    let record = tampered_records
        .iter_mut()
        .find(|record| record.id == reference.id)
        .expect("built-in count evidence");
    record.datum = EvidenceDatum::Count(95);
    let tampered = BenchmarkCorpus {
        retained_evidence: &tampered_records,
        ..benchmark_corpus()
    };
    assert!(matches!(
        resolve_evidence(&tampered, reference),
        Err(EvidenceError::Tampered { .. })
    ));

    let original = retained_evidence()[3].clone();
    let mut alias = original.clone();
    alias.id = "evidence:alias-with-same-authoritative-content:v1";
    let duplicate_digest_records = [original.clone(), alias];
    let duplicate_digest_corpus = BenchmarkCorpus {
        retained_evidence: &duplicate_digest_records,
        ..benchmark_corpus()
    };
    assert!(matches!(
        resolve_evidence(
            &duplicate_digest_corpus,
            EvidenceRef::exact(original.id, original.semantic_digest())
        ),
        Err(EvidenceError::DuplicateDigest { .. })
    ));
    assert!(
        audit_corpus(&duplicate_digest_corpus)
            .gaps
            .iter()
            .any(|gap| gap.contains("duplicate retained evidence digest"))
    );
}

#[test]
fn typed_rate_refuses_zero_denominator_and_invalid_count_domain() {
    const EVALUATOR: &str = "count-rate-v1";
    const PROPOSAL: &str = "test-proposal-v1";
    let numerator = EvidenceRecord::count(
        "count:numerator",
        "subject",
        "events",
        EVALUATOR,
        PROPOSAL,
        EvidenceRole::RateNumerator,
        1,
    );
    let zero = EvidenceRecord::count(
        "count:denominator",
        "subject",
        "attempts",
        EVALUATOR,
        PROPOSAL,
        EvidenceRole::RateDenominator,
        0,
    );
    let request = MetricRequest {
        schema: MetricSchema::GenericRate,
        evaluator_semantics: EVALUATOR,
        proposal_semantics: PROPOSAL,
        definition: MetricDefinition::Rate {
            numerator: role(&numerator, "subject", "events", "count"),
            denominator: role(&zero, "subject", "attempts", "count"),
        },
    };
    let records = [numerator.clone(), zero];
    let corpus = BenchmarkCorpus {
        retained_evidence: &records,
        ..benchmark_corpus()
    };
    assert_eq!(
        reconstruct_metric(&corpus, request),
        Err(MetricError::ZeroDenominator)
    );

    let denominator = EvidenceRecord::count(
        "count:denominator",
        "subject",
        "attempts",
        EVALUATOR,
        PROPOSAL,
        EvidenceRole::RateDenominator,
        1,
    );
    let excess = EvidenceRecord::count(
        "count:numerator",
        "subject",
        "events",
        EVALUATOR,
        PROPOSAL,
        EvidenceRole::RateNumerator,
        2,
    );
    let records = [excess.clone(), denominator.clone()];
    let corpus = BenchmarkCorpus {
        retained_evidence: &records,
        ..benchmark_corpus()
    };
    let excess_request = MetricRequest {
        definition: MetricDefinition::Rate {
            numerator: role(&excess, "subject", "events", "count"),
            denominator: role(&denominator, "subject", "attempts", "count"),
        },
        ..request
    };
    assert_eq!(
        reconstruct_metric(&corpus, excess_request),
        Err(MetricError::NumeratorExceedsDenominator)
    );

    let swap_records = [numerator.clone(), denominator.clone()];
    let swap_corpus = BenchmarkCorpus {
        retained_evidence: &swap_records,
        ..benchmark_corpus()
    };
    let swapped = MetricRequest {
        definition: MetricDefinition::Rate {
            numerator: role(&denominator, "subject", "attempts", "count"),
            denominator: role(&numerator, "subject", "events", "count"),
        },
        ..request
    };
    assert_eq!(
        reconstruct_metric(&swap_corpus, swapped),
        Err(MetricError::RoleMismatch {
            expected: EvidenceRole::RateNumerator,
            actual: EvidenceRole::RateDenominator,
        })
    );

    let large_numerator = EvidenceRecord::count(
        "count:numerator",
        "subject",
        "events",
        EVALUATOR,
        PROPOSAL,
        EvidenceRole::RateNumerator,
        1_u64 << 53,
    );
    let large_denominator = EvidenceRecord::count(
        "count:denominator",
        "subject",
        "attempts",
        EVALUATOR,
        PROPOSAL,
        EvidenceRole::RateDenominator,
        (1_u64 << 53) + 1,
    );
    let large_records = [large_numerator.clone(), large_denominator.clone()];
    let large_corpus = BenchmarkCorpus {
        retained_evidence: &large_records,
        ..benchmark_corpus()
    };
    let large_request = MetricRequest {
        definition: MetricDefinition::Rate {
            numerator: role(&large_numerator, "subject", "events", "count"),
            denominator: role(&large_denominator, "subject", "attempts", "count"),
        },
        ..request
    };
    assert_eq!(
        reconstruct_metric(&large_corpus, large_request),
        Err(MetricError::InexactCount {
            role: "denominator"
        })
    );
}

#[test]
fn governance_rows_bind_typed_schemas_and_missing_evidence_is_unavailable() {
    let proposals = instrumented_proposals();
    assert_eq!(proposals.len(), 8);
    let ids: Vec<_> = proposals.iter().map(|proposal| proposal.proposal).collect();
    for id in ["8", "1", "2", "D", "F", "9", "10", "A"] {
        assert!(ids.contains(&id), "proposal {id} not instrumented");
    }
    for proposal in proposals {
        assert!(!proposal.dataset.name().is_empty());
        assert!(!proposal.evaluator.formulas().is_empty());
        assert!(!proposal.evaluator.required_roles().is_empty());
        assert!(!proposal.kill_metric.is_empty());
        assert!(!proposal.evaluator_semantics.is_empty());
        assert!(!proposal.proposal_semantics.is_empty());
        let ProposalEvaluation::Unavailable { missing } =
            evaluate_proposal(&benchmark_corpus(), proposal)
        else {
            panic!(
                "built-in proposal {} must lack real retained measurements",
                proposal.proposal
            );
        };
        assert!(!missing.is_empty());
        assert!(
            missing
                .iter()
                .all(|entry| entry.proposal == proposal.proposal)
        );
    }
    let ProposalEvaluation::Unavailable { missing } =
        evaluate_proposal(&benchmark_corpus(), &proposals[0])
    else {
        unreachable!("proposal 8 has no retained planner evidence");
    };
    assert_eq!(missing.len(), query_set().len() * 2);
    assert_eq!(missing[0].subject_id, query_set()[0].id);
    assert_eq!(missing[0].role, EvidenceRole::PlannerBaselineCost);
    assert_eq!(missing[1].role, EvidenceRole::PlannerCandidateCost);

    let mut mismatched = proposals[0];
    mismatched.dataset = DatasetKind::DesignTasks;
    assert!(matches!(
        evaluate_proposal(&benchmark_corpus(), &mismatched),
        ProposalEvaluation::Refused { diagnostics }
            if diagnostics.contains(&ProposalRefusal::SchemaMismatch {
                proposal: "8",
                field: "dataset",
            })
    ));
    let no_schemas = audit_corpus(&BenchmarkCorpus {
        instrumented_proposals: &[],
        ..benchmark_corpus()
    });
    assert!(!no_schemas.ok());
    assert!(
        no_schemas
            .gaps
            .iter()
            .any(|gap| gap.contains("proposal evaluator schema count"))
    );
    let mut duplicate_queries = query_set().to_vec();
    let duplicate_id = duplicate_queries[0].id;
    duplicate_queries[1].id = duplicate_id;
    let duplicate_corpus = BenchmarkCorpus {
        query_set: &duplicate_queries,
        ..benchmark_corpus()
    };
    assert!(matches!(
        evaluate_proposal(&duplicate_corpus, &proposals[0]),
        ProposalEvaluation::Refused { diagnostics }
            if matches!(
                diagnostics.as_slice(),
                [ProposalRefusal::DuplicateSubject { proposal: "8", .. }]
            )
    ));
}

#[test]
fn complete_typed_proposal_evidence_reconstructs_available_metrics() {
    let proposal = &instrumented_proposals()[0];
    let ids = [
        ("p8:q1:baseline", "p8:q1:candidate"),
        ("p8:q2:baseline", "p8:q2:candidate"),
        ("p8:q3:baseline", "p8:q3:candidate"),
    ];
    let mut records = Vec::new();
    for (query, (baseline_id, candidate_id)) in query_set().iter().zip(ids) {
        records.push(EvidenceRecord::scalar(
            baseline_id,
            query.id,
            "planner-baseline-cost",
            "work-unit",
            proposal.evaluator_semantics,
            proposal.proposal_semantics,
            EvidenceRole::PlannerBaselineCost,
            100.0,
        ));
        records.push(EvidenceRecord::scalar(
            candidate_id,
            query.id,
            "planner-candidate-cost",
            "work-unit",
            proposal.evaluator_semantics,
            proposal.proposal_semantics,
            EvidenceRole::PlannerCandidateCost,
            50.0,
        ));
    }
    let manifest: Vec<_> = records
        .iter()
        .map(|record| ProposalEvidenceRef {
            proposal: proposal.proposal,
            subject_id: record.subject_id,
            role: record.role,
            reference: EvidenceRef::exact(record.id, record.semantic_digest()),
        })
        .collect();
    let unreferenced = BenchmarkCorpus {
        retained_evidence: &records,
        ..benchmark_corpus()
    };
    assert!(matches!(
        evaluate_proposal(&unreferenced, proposal),
        ProposalEvaluation::Unavailable { .. }
    ));
    let corpus = BenchmarkCorpus {
        retained_evidence: &records,
        proposal_evidence: &manifest,
        ..benchmark_corpus()
    };
    let ProposalEvaluation::Available { metrics } = evaluate_proposal(&corpus, proposal) else {
        panic!("complete typed evidence must reconstruct");
    };
    assert_eq!(metrics.len(), query_set().len());
    assert!(metrics.iter().all(|metric| metric.value == 2.0));
    assert_eq!(
        ProposalEvaluation::Available {
            metrics: metrics.clone()
        },
        evaluate_proposal(&corpus, proposal)
    );

    let mut tampered_manifest = manifest.clone();
    tampered_manifest[0].reference.expected_digest = ContentHash([9; 32]);
    let tampered = BenchmarkCorpus {
        retained_evidence: &records,
        proposal_evidence: &tampered_manifest,
        ..benchmark_corpus()
    };
    assert!(matches!(
        evaluate_proposal(&tampered, proposal),
        ProposalEvaluation::Refused { diagnostics }
            if matches!(
                diagnostics.as_slice(),
                [ProposalRefusal::Evidence {
                    error: EvidenceError::Tampered { .. },
                    ..
                }]
            )
    ));
    let mut wrong_units_records = records.clone();
    wrong_units_records[0].units = "seconds";
    let mut resigned_manifest = manifest.clone();
    resigned_manifest[0].reference = EvidenceRef::exact(
        wrong_units_records[0].id,
        wrong_units_records[0].semantic_digest(),
    );
    let wrong_units = BenchmarkCorpus {
        retained_evidence: &wrong_units_records,
        proposal_evidence: &resigned_manifest,
        ..benchmark_corpus()
    };
    assert!(matches!(
        evaluate_proposal(&wrong_units, proposal),
        ProposalEvaluation::Refused { diagnostics }
            if matches!(
                diagnostics.as_slice(),
                [ProposalRefusal::Evidence {
                    error: EvidenceError::ContextMismatch { field: "units", .. },
                    ..
                }]
            )
    ));
}

#[test]
fn special_proposal_formulas_reconstruct_from_independent_manifests() {
    let merge = instrumented_proposals()
        .iter()
        .find(|proposal| proposal.evaluator == ProposalEvaluator::Merge10)
        .expect("merge schema");
    let merge_subject = merge_trials()[0].id;
    let merge_records = [
        EvidenceRecord::count(
            "merge10:conflicts",
            merge_subject,
            EvidenceRole::MergeConflicts.name(),
            merge.evaluator_semantics,
            merge.proposal_semantics,
            EvidenceRole::MergeConflicts,
            1,
        ),
        EvidenceRecord::count(
            "merge10:escalations",
            merge_subject,
            EvidenceRole::MergeEscalations.name(),
            merge.evaluator_semantics,
            merge.proposal_semantics,
            EvidenceRole::MergeEscalations,
            1,
        ),
        EvidenceRecord::count(
            "merge10:refusals",
            merge_subject,
            EvidenceRole::MergeRefusals.name(),
            merge.evaluator_semantics,
            merge.proposal_semantics,
            EvidenceRole::MergeRefusals,
            0,
        ),
        EvidenceRecord::count(
            "merge10:type-conflicts",
            merge_subject,
            EvidenceRole::MergeTypeConflicts.name(),
            merge.evaluator_semantics,
            merge.proposal_semantics,
            EvidenceRole::MergeTypeConflicts,
            0,
        ),
        EvidenceRecord::count(
            "merge10:attempts",
            merge_subject,
            EvidenceRole::MergeAttempts.name(),
            merge.evaluator_semantics,
            merge.proposal_semantics,
            EvidenceRole::MergeAttempts,
            10,
        ),
    ];
    let merge_manifest = proposal_manifest(merge, &merge_records);
    assert_proposal_values(
        &BenchmarkCorpus {
            merge_trials: &merge_trials()[..1],
            retained_evidence: &merge_records,
            proposal_evidence: &merge_manifest,
            ..benchmark_corpus()
        },
        merge,
        &[0.2],
    );

    let coverage = instrumented_proposals()
        .iter()
        .find(|proposal| proposal.evaluator == ProposalEvaluator::CoverageA)
        .expect("coverage schema");
    let coverage_subject = query_set()[0].id;
    let coverage_records = [
        EvidenceRecord::scalar(
            "coverage-a:covered",
            coverage_subject,
            EvidenceRole::CoveredQueryVolume.name(),
            "query-volume",
            coverage.evaluator_semantics,
            coverage.proposal_semantics,
            EvidenceRole::CoveredQueryVolume,
            2.0,
        ),
        EvidenceRecord::scalar(
            "coverage-a:total",
            coverage_subject,
            EvidenceRole::TotalQueryVolume.name(),
            "query-volume",
            coverage.evaluator_semantics,
            coverage.proposal_semantics,
            EvidenceRole::TotalQueryVolume,
            10.0,
        ),
    ];
    let coverage_manifest = proposal_manifest(coverage, &coverage_records);
    assert_proposal_values(
        &BenchmarkCorpus {
            query_set: &query_set()[..1],
            retained_evidence: &coverage_records,
            proposal_evidence: &coverage_manifest,
            ..benchmark_corpus()
        },
        coverage,
        &[0.2],
    );

    let guard = instrumented_proposals()
        .iter()
        .find(|proposal| proposal.evaluator == ProposalEvaluator::GuardD)
        .expect("guard schema");
    let guard_subject = design_tasks()[0].id;
    let guard_records = [
        EvidenceRecord::count(
            "guard-d:endpoint-catches",
            guard_subject,
            EvidenceRole::GuardEndpointCatches.name(),
            guard.evaluator_semantics,
            guard.proposal_semantics,
            EvidenceRole::GuardEndpointCatches,
            8,
        ),
        EvidenceRecord::count(
            "guard-d:endpoint-trials",
            guard_subject,
            EvidenceRole::GuardEndpointTrials.name(),
            guard.evaluator_semantics,
            guard.proposal_semantics,
            EvidenceRole::GuardEndpointTrials,
            10,
        ),
        EvidenceRecord::count(
            "guard-d:random-catches",
            guard_subject,
            EvidenceRole::GuardRandomCatches.name(),
            guard.evaluator_semantics,
            guard.proposal_semantics,
            EvidenceRole::GuardRandomCatches,
            3,
        ),
        EvidenceRecord::count(
            "guard-d:random-trials",
            guard_subject,
            EvidenceRole::GuardRandomTrials.name(),
            guard.evaluator_semantics,
            guard.proposal_semantics,
            EvidenceRole::GuardRandomTrials,
            10,
        ),
    ];
    let guard_manifest = proposal_manifest(guard, &guard_records);
    assert_proposal_values(
        &BenchmarkCorpus {
            design_tasks: &design_tasks()[..1],
            retained_evidence: &guard_records,
            proposal_evidence: &guard_manifest,
            ..benchmark_corpus()
        },
        guard,
        &[0.8, 0.3],
    );

    let speculation = instrumented_proposals()
        .iter()
        .find(|proposal| proposal.evaluator == ProposalEvaluator::Speculation9)
        .expect("speculation schema");
    let speculation_subject = mms_battery()[0].id;
    let speculation_records = [
        EvidenceRecord::count(
            "speculation9:accepts",
            speculation_subject,
            EvidenceRole::SpeculationAccepts.name(),
            speculation.evaluator_semantics,
            speculation.proposal_semantics,
            EvidenceRole::SpeculationAccepts,
            4,
        ),
        EvidenceRecord::count(
            "speculation9:attempts",
            speculation_subject,
            EvidenceRole::SpeculationAttempts.name(),
            speculation.evaluator_semantics,
            speculation.proposal_semantics,
            EvidenceRole::SpeculationAttempts,
            10,
        ),
        EvidenceRecord::scalar(
            "speculation9:cold-cost",
            speculation_subject,
            EvidenceRole::ColdStartCost.name(),
            "work-unit",
            speculation.evaluator_semantics,
            speculation.proposal_semantics,
            EvidenceRole::ColdStartCost,
            30.0,
        ),
        EvidenceRecord::scalar(
            "speculation9:warm-cost",
            speculation_subject,
            EvidenceRole::WarmStartCost.name(),
            "work-unit",
            speculation.evaluator_semantics,
            speculation.proposal_semantics,
            EvidenceRole::WarmStartCost,
            10.0,
        ),
    ];
    let speculation_manifest = proposal_manifest(speculation, &speculation_records);
    assert_proposal_values(
        &BenchmarkCorpus {
            mms_battery: &mms_battery()[..1],
            retained_evidence: &speculation_records,
            proposal_evidence: &speculation_manifest,
            ..benchmark_corpus()
        },
        speculation,
        &[0.4, 3.0],
    );
}

#[test]
fn corpus_identity_is_deterministic_length_framed_and_schema_versioned() {
    assert_eq!(corpus_digest(), corpus_digest());
    assert_ne!(corpus_digest(), ContentHash([0; 32]));

    let mut left_rows = query_set().to_vec();
    left_rows[0].id = "ab";
    left_rows[0].qoi = "c";
    let left = BenchmarkCorpus {
        query_set: &left_rows,
        ..benchmark_corpus()
    };
    let mut right_rows = query_set().to_vec();
    right_rows[0].id = "a";
    right_rows[0].qoi = "bc";
    let right = BenchmarkCorpus {
        query_set: &right_rows,
        ..benchmark_corpus()
    };
    assert_ne!(corpus_digest_for(&left), corpus_digest_for(&right));

    assert_digest_changed(&BenchmarkCorpus {
        version: BENCHMARK_VERSION + 1,
        ..benchmark_corpus()
    });
    assert_digest_changed(&BenchmarkCorpus {
        identity_schema_version: CORPUS_IDENTITY_SCHEMA_VERSION + 1,
        ..benchmark_corpus()
    });
}

#[test]
#[allow(clippy::too_many_lines)]
fn every_corpus_semantic_field_is_mutation_sensitive() {
    macro_rules! mutate_row {
        ($corpus_field:ident, $source:expr, $row_field:ident, $value:expr) => {{
            let mut rows = ($source).to_vec();
            rows[0].$row_field = $value;
            let changed = BenchmarkCorpus {
                $corpus_field: &rows,
                ..benchmark_corpus()
            };
            assert_digest_changed(&changed);
        }};
    }

    mutate_row!(query_set, query_set(), id, "changed-query");
    mutate_row!(query_set, query_set(), qoi, "changed QoI");
    mutate_row!(query_set, query_set(), units, "changed-unit");
    mutate_row!(query_set, query_set(), tolerance, 0.25);
    mutate_row!(query_set, query_set(), reference_answer, 358.25);
    mutate_row!(query_set, query_set(), reference_cost, 1001.0);
    mutate_row!(
        query_set,
        query_set(),
        reference_cost_units,
        "changed-cost-unit"
    );
    mutate_row!(
        query_set,
        query_set(),
        reference_evidence,
        EvidenceRef::exact("changed-evidence", ContentHash([1; 32]))
    );
    mutate_row!(
        query_set,
        query_set(),
        reference_evidence,
        EvidenceRef {
            expected_digest: changed_digest(query_set()[0].reference_evidence.expected_digest),
            ..query_set()[0].reference_evidence
        }
    );
    mutate_row!(
        query_set,
        query_set(),
        reference_evaluator_semantics,
        "changed-evaluator"
    );
    mutate_row!(
        query_set,
        query_set(),
        reference_proposal_semantics,
        "changed-proposal"
    );

    mutate_row!(design_tasks, design_tasks(), id, "changed-task");
    mutate_row!(design_tasks, design_tasks(), dimension, 5);
    mutate_row!(design_tasks, design_tasks(), optimum, 350.5);

    mutate_row!(edit_traces, edit_traces(), id, "changed-trace");
    mutate_row!(edit_traces, edit_traces(), total_ops, 121);
    mutate_row!(edit_traces, edit_traces(), correct_skips, 95);
    mutate_row!(
        edit_traces,
        edit_traces(),
        correct_skips_evidence,
        EvidenceRef::exact("changed-skip-evidence", ContentHash([1; 32]))
    );
    mutate_row!(
        edit_traces,
        edit_traces(),
        total_ops_evidence,
        EvidenceRef::exact("changed-total-evidence", ContentHash([1; 32]))
    );
    mutate_row!(
        edit_traces,
        edit_traces(),
        evaluator_semantics,
        "changed-edit-evaluator"
    );
    mutate_row!(
        edit_traces,
        edit_traces(),
        proposal_semantics,
        "changed-edit-proposal"
    );

    mutate_row!(mms_battery, mms_battery(), id, "changed-mms");
    mutate_row!(mms_battery, mms_battery(), exact_center, 0.75);

    mutate_row!(merge_trials, merge_trials(), id, "changed-merge");
    mutate_row!(merge_trials, merge_trials(), total_merges, 41);
    mutate_row!(
        merge_trials,
        merge_trials(),
        candidate_remainder_conflicts,
        7
    );
    mutate_row!(
        merge_trials,
        merge_trials(),
        conflict_count_evidence,
        EvidenceRef::exact("changed-conflict-evidence", ContentHash([1; 32]))
    );
    mutate_row!(
        merge_trials,
        merge_trials(),
        total_merges_evidence,
        EvidenceRef::exact("changed-merge-total-evidence", ContentHash([1; 32]))
    );
    mutate_row!(
        merge_trials,
        merge_trials(),
        evaluator_semantics,
        "changed-merge-evaluator"
    );
    mutate_row!(
        merge_trials,
        merge_trials(),
        proposal_semantics,
        "changed-merge-proposal"
    );

    mutate_row!(retained_evidence, retained_evidence(), id, "changed-record");
    mutate_row!(
        retained_evidence,
        retained_evidence(),
        subject_id,
        "changed-subject"
    );
    mutate_row!(
        retained_evidence,
        retained_evidence(),
        quantity,
        "changed-quantity"
    );
    mutate_row!(
        retained_evidence,
        retained_evidence(),
        units,
        "changed-unit"
    );
    mutate_row!(
        retained_evidence,
        retained_evidence(),
        role,
        EvidenceRole::BaselineCost
    );
    mutate_row!(
        retained_evidence,
        retained_evidence(),
        evaluator_semantics,
        "changed-evidence-evaluator"
    );
    mutate_row!(
        retained_evidence,
        retained_evidence(),
        proposal_semantics,
        "changed-evidence-proposal"
    );
    macro_rules! mutate_reference_field {
        ($field:ident, $value:expr) => {{
            let mut rows = retained_evidence().to_vec();
            let EvidenceDatum::QueryReference { $field, .. } = &mut rows[0].datum else {
                panic!("first retained record must be query evidence");
            };
            *$field = $value;
            assert_digest_changed(&BenchmarkCorpus {
                retained_evidence: &rows,
                ..benchmark_corpus()
            });
        }};
    }
    mutate_reference_field!(answer, 358.25);
    mutate_reference_field!(tolerance, 0.25);
    mutate_reference_field!(reference_cost, 1001.0);
    mutate_reference_field!(reference_cost_units, "changed-cost-unit");
    mutate_reference_field!(
        color,
        Color::Estimated {
            estimator: "changed-estimator-v1".to_string(),
            dispersion: 1.0,
        }
    );
    mutate_reference_field!(
        admission_receipt,
        Some(AdmissionReceipt::from_parts(
            ContentHash([3; 32]),
            1,
            COLOR_ALGEBRA_VERSION,
            ContentHash([4; 32]),
        ))
    );
    mutate_row!(
        retained_evidence,
        retained_evidence(),
        datum,
        EvidenceDatum::Count(358)
    );
    mutate_row!(
        instrumented_proposals,
        instrumented_proposals(),
        proposal,
        "changed-id"
    );
    mutate_row!(
        instrumented_proposals,
        instrumented_proposals(),
        dataset,
        DatasetKind::DesignTasks
    );
    mutate_row!(
        instrumented_proposals,
        instrumented_proposals(),
        evaluator,
        ProposalEvaluator::Adjoint1
    );
    mutate_row!(
        instrumented_proposals,
        instrumented_proposals(),
        kill_metric,
        "changed-kill-metric"
    );
    mutate_row!(
        instrumented_proposals,
        instrumented_proposals(),
        evaluator_semantics,
        "changed-governance-evaluator"
    );
    mutate_row!(
        instrumented_proposals,
        instrumented_proposals(),
        proposal_semantics,
        "changed-governance-proposal"
    );

    let base_manifest = ProposalEvidenceRef {
        proposal: "8",
        subject_id: "cht-q1",
        role: EvidenceRole::PlannerBaselineCost,
        reference: EvidenceRef::exact("manifest-evidence", ContentHash([5; 32])),
    };
    let base_entries = [base_manifest];
    let base_manifest_digest = corpus_digest_for(&BenchmarkCorpus {
        proposal_evidence: &base_entries,
        ..benchmark_corpus()
    });
    for changed in [
        ProposalEvidenceRef {
            proposal: "1",
            ..base_manifest
        },
        ProposalEvidenceRef {
            subject_id: "cht-q2",
            ..base_manifest
        },
        ProposalEvidenceRef {
            role: EvidenceRole::PlannerCandidateCost,
            ..base_manifest
        },
        ProposalEvidenceRef {
            reference: EvidenceRef::exact("changed-manifest-evidence", ContentHash([5; 32])),
            ..base_manifest
        },
        ProposalEvidenceRef {
            reference: EvidenceRef::exact("manifest-evidence", ContentHash([6; 32])),
            ..base_manifest
        },
    ] {
        let entries = [changed];
        assert_ne!(
            corpus_digest_for(&BenchmarkCorpus {
                proposal_evidence: &entries,
                ..benchmark_corpus()
            }),
            base_manifest_digest
        );
    }

    let mut reordered_queries = query_set().to_vec();
    reordered_queries.swap(0, 1);
    assert_digest_changed(&BenchmarkCorpus {
        query_set: &reordered_queries,
        ..benchmark_corpus()
    });
    assert_digest_changed(&BenchmarkCorpus {
        query_set: &query_set()[1..],
        ..benchmark_corpus()
    });
    assert_digest_changed(&BenchmarkCorpus {
        design_tasks: &design_tasks()[1..],
        ..benchmark_corpus()
    });
    assert_digest_changed(&BenchmarkCorpus {
        edit_traces: &edit_traces()[1..],
        ..benchmark_corpus()
    });
    assert_digest_changed(&BenchmarkCorpus {
        mms_battery: &mms_battery()[1..],
        ..benchmark_corpus()
    });
    assert_digest_changed(&BenchmarkCorpus {
        merge_trials: &merge_trials()[1..],
        ..benchmark_corpus()
    });
    assert_digest_changed(&BenchmarkCorpus {
        retained_evidence: &retained_evidence()[1..],
        ..benchmark_corpus()
    });
    assert_digest_changed(&BenchmarkCorpus {
        instrumented_proposals: &instrumented_proposals()[1..],
        ..benchmark_corpus()
    });
}

#[test]
fn audit_rejects_duplicate_ids_nonfinite_domains_and_invalid_counts() {
    let mut duplicate_queries = query_set().to_vec();
    let duplicate_id = duplicate_queries[0].id;
    duplicate_queries[1].id = duplicate_id;
    let duplicate = audit_corpus(&BenchmarkCorpus {
        query_set: &duplicate_queries,
        ..benchmark_corpus()
    });
    assert!(
        duplicate
            .gaps
            .iter()
            .any(|gap| gap.contains("duplicate corpus row id"))
    );

    let mut invalid_queries = query_set().to_vec();
    invalid_queries[0].tolerance = f64::NAN;
    invalid_queries[1].reference_answer = f64::INFINITY;
    invalid_queries[2].reference_cost = 0.0;
    let invalid_query_audit = audit_corpus(&BenchmarkCorpus {
        query_set: &invalid_queries,
        ..benchmark_corpus()
    });
    assert!(
        invalid_query_audit
            .gaps
            .iter()
            .any(|gap| gap.contains("invalid tolerance"))
    );
    assert!(
        invalid_query_audit
            .gaps
            .iter()
            .any(|gap| gap.contains("non-finite reference answer"))
    );
    assert!(
        invalid_query_audit
            .gaps
            .iter()
            .any(|gap| gap.contains("invalid reference cost"))
    );

    let mut invalid_tasks = design_tasks().to_vec();
    invalid_tasks[0].dimension = 0;
    invalid_tasks[1].optimum = f64::NAN;
    let task_audit = audit_corpus(&BenchmarkCorpus {
        design_tasks: &invalid_tasks,
        ..benchmark_corpus()
    });
    assert!(
        task_audit
            .gaps
            .iter()
            .any(|gap| gap.contains("zero dimension"))
    );
    assert!(
        task_audit
            .gaps
            .iter()
            .any(|gap| gap.contains("non-finite optimum"))
    );

    let mut invalid_edits = edit_traces().to_vec();
    invalid_edits[0].total_ops = 0;
    invalid_edits[0].correct_skips = 1;
    let edit_audit = audit_corpus(&BenchmarkCorpus {
        edit_traces: &invalid_edits,
        ..benchmark_corpus()
    });
    assert!(
        edit_audit
            .gaps
            .iter()
            .any(|gap| gap.contains("zero total ops"))
    );
    assert!(
        edit_audit
            .gaps
            .iter()
            .any(|gap| gap.contains("skips exceed total ops"))
    );

    let mut invalid_merges = merge_trials().to_vec();
    invalid_merges[0].total_merges = 0;
    invalid_merges[0].candidate_remainder_conflicts = 1;
    let merge_audit = audit_corpus(&BenchmarkCorpus {
        merge_trials: &invalid_merges,
        ..benchmark_corpus()
    });
    assert!(
        merge_audit
            .gaps
            .iter()
            .any(|gap| gap.contains("zero total merges"))
    );
    assert!(
        merge_audit
            .gaps
            .iter()
            .any(|gap| gap.contains("conflicts exceed total merges"))
    );
}

#[test]
fn audit_rejects_duplicate_or_tampered_evidence_and_governance_semantics() {
    let mut duplicate_evidence = retained_evidence().to_vec();
    duplicate_evidence.push(retained_evidence()[0].clone());
    let duplicate = audit_corpus(&BenchmarkCorpus {
        retained_evidence: &duplicate_evidence,
        ..benchmark_corpus()
    });
    assert!(
        duplicate
            .gaps
            .iter()
            .any(|gap| gap.contains("duplicate retained evidence id"))
    );

    let mut tampered_evidence = retained_evidence().to_vec();
    tampered_evidence[0].quantity = "tampered";
    let tampered = audit_corpus(&BenchmarkCorpus {
        retained_evidence: &tampered_evidence,
        ..benchmark_corpus()
    });
    assert!(
        tampered
            .gaps
            .iter()
            .any(|gap| gap.contains("mismatches its expected content identity"))
    );

    let nonfinite_record = EvidenceRecord::query_reference(
        retained_evidence()[0].id,
        query_set()[0].id,
        query_set()[0].qoi,
        query_set()[0].units,
        query_set()[0].reference_evaluator_semantics,
        query_set()[0].reference_proposal_semantics,
        f64::NAN,
        query_set()[0].tolerance,
        query_set()[0].reference_cost,
        query_set()[0].reference_cost_units,
        Color::Estimated {
            estimator: "nonfinite-query-v1".to_string(),
            dispersion: f64::INFINITY,
        },
        None,
    );
    let nonfinite_records = [nonfinite_record];
    let nonfinite = audit_corpus(&BenchmarkCorpus {
        retained_evidence: &nonfinite_records,
        ..benchmark_corpus()
    });
    assert!(
        nonfinite
            .gaps
            .iter()
            .any(|gap| gap.contains("non-finite scalar"))
    );

    let mut bad_governance: Vec<InstrumentedProposal> = instrumented_proposals().to_vec();
    bad_governance[0].evaluator_semantics = "";
    let duplicate_proposal = bad_governance[0].proposal;
    bad_governance[1].proposal = duplicate_proposal;
    let governance = audit_corpus(&BenchmarkCorpus {
        instrumented_proposals: &bad_governance,
        ..benchmark_corpus()
    });
    assert!(
        governance
            .gaps
            .iter()
            .any(|gap| gap.contains("empty evaluator semantics"))
    );
    assert!(
        governance
            .gaps
            .iter()
            .any(|gap| gap.contains("duplicate instrumented proposal"))
    );
}

#[test]
fn built_in_audit_fails_closed_until_real_proposal_evidence_exists() {
    let result = audit();
    assert!(!result.ok());
    assert_eq!(result.version, BENCHMARK_VERSION);
    assert_eq!(result.instrumented, 0);
    assert!(
        result
            .gaps
            .iter()
            .any(|gap| gap.contains("proposal 8 unavailable"))
    );
    assert!(
        result
            .gaps
            .iter()
            .any(|gap| gap.contains("proposal 10 unavailable"))
    );
}
