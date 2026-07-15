//! G0/G3/G4/G5 conformance tests for RE.Z1 hybrid-time and Zeno semantics.

#![allow(clippy::wildcard_imports)]
#![allow(
    clippy::too_many_lines,
    reason = "each long RE.Z1 case keeps one hybrid semantics and refusal narrative auditable"
)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_time::hybrid::*;

fn bytes(seed: u8) -> [u8; 32] {
    [seed; 32]
}

fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x5A45_4E4F,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn mode(seed: u8) -> HybridModeSpecV1 {
    HybridModeSpecV1 {
        mode: HybridModeIdV1::from_bytes(bytes(seed)),
        dynamics: ContinuousDynamicsIdV1::from_bytes(bytes(seed.wrapping_add(1))),
        state_dimension: 2,
        class: ContinuousDynamicsClassV1::DeterministicOde,
    }
}

fn event(
    seed: u8,
    source: HybridModeIdV1,
    target: HybridModeIdV1,
    law: InteractionLawV1,
    dwell: DwellSemanticsV1,
) -> HybridEventSpecV1 {
    HybridEventSpecV1 {
        event: HybridEventIdV1::from_bytes(bytes(seed)),
        source_mode: source,
        guard: HybridGuardIdV1::from_bytes(bytes(seed.wrapping_add(1))),
        orientation: GuardOrientationV1::NegativeToPositive,
        crossing: CrossingSemanticsV1::Transverse {
            witness: HybridWitnessIdV1::from_bytes(bytes(seed.wrapping_add(2))),
        },
        reset: ResetSemanticsV1::Deterministic {
            relation: ResetRelationIdV1::from_bytes(bytes(seed.wrapping_add(3))),
            target,
        },
        law,
        simultaneity: EventSimultaneityV1::Exclusive {
            witness: HybridWitnessIdV1::from_bytes(bytes(seed.wrapping_add(4))),
        },
        dwell,
    }
}

fn relay_ir() -> ZenoProblemIrV1 {
    let left = mode(10);
    let right = mode(20);
    let left_to_right = event(
        30,
        left.mode,
        right.mode,
        InteractionLawV1::Relay {
            law: InteractionLawIdV1::from_bytes(bytes(31)),
        },
        DwellSemanticsV1::ZeroAllowed,
    );
    let right_to_left = event(
        40,
        right.mode,
        left.mode,
        InteractionLawV1::Relay {
            law: InteractionLawIdV1::from_bytes(bytes(41)),
        },
        DwellSemanticsV1::Unknown {
            no_claim: HybridNoClaimIdV1::from_bytes(bytes(42)),
        },
    );
    ZenoProblemIrV1::new(
        HybridModelIdV1::from_bytes(bytes(1)),
        HybridModelVersionIdV1::from_bytes(bytes(2)),
        HybridModelLineageV1::Original,
        HybridFrameIdV1::from_bytes(bytes(3)),
        HybridUnitSystemIdV1::from_bytes(bytes(4)),
        HybridStateSetIdV1::from_bytes(bytes(5)),
        HybridTimeScaleV1 {
            unit: HybridTimeUnitIdV1::from_bytes(bytes(6)),
            seconds_per_unit: 1.0,
        },
        vec![right, left],
        vec![right_to_left, left_to_right],
        EventLanguageSpecV1 {
            language: EventLanguageIdV1::from_bytes(bytes(7)),
            semantics: EventLanguageSemanticsV1::OmegaLanguage,
        },
        SimultaneousEventPolicyV1::NoSimultaneousEvents {
            witness: HybridWitnessIdV1::from_bytes(bytes(8)),
        },
        HybridTimeDomainV1 {
            start: 0.0,
            end: HybridTimeEndV1::Finite(10.0),
            event_cap: Some(10_000),
        },
        CompactnessSemanticsV1::Compact {
            witness: HybridWitnessIdV1::from_bytes(bytes(9)),
        },
        AccumulationCandidateV1::Window {
            earliest: 0.9,
            latest: 1.1,
            states: HybridStateSetIdV1::from_bytes(bytes(50)),
            trace: HybridEventTraceIdV1::from_bytes(bytes(51)),
        },
        ContinuationCategoryV1::SetValued {
            rule: ContinuationRuleIdV1::from_bytes(bytes(52)),
        },
        HybridAnalysisBudgetV1 {
            max_event_word_len: 100_000,
            max_transitions: 1_000_000,
            max_wall_seconds: 5.0,
        },
    )
}

fn separated_ir() -> ZenoProblemIrV1 {
    let before = mode(60);
    let after = mode(70);
    let transition = event(
        80,
        before.mode,
        after.mode,
        InteractionLawV1::None {
            justification: HybridNoClaimIdV1::from_bytes(bytes(81)),
        },
        DwellSemanticsV1::PositiveLowerBound {
            value: 0.25,
            witness: HybridWitnessIdV1::from_bytes(bytes(82)),
        },
    );
    ZenoProblemIrV1::new(
        HybridModelIdV1::from_bytes(bytes(53)),
        HybridModelVersionIdV1::from_bytes(bytes(54)),
        HybridModelLineageV1::Original,
        HybridFrameIdV1::from_bytes(bytes(55)),
        HybridUnitSystemIdV1::from_bytes(bytes(56)),
        HybridStateSetIdV1::from_bytes(bytes(57)),
        HybridTimeScaleV1 {
            unit: HybridTimeUnitIdV1::from_bytes(bytes(58)),
            seconds_per_unit: 1.0,
        },
        vec![before, after],
        vec![transition],
        EventLanguageSpecV1 {
            language: EventLanguageIdV1::from_bytes(bytes(59)),
            semantics: EventLanguageSemanticsV1::FiniteWords {
                max_events_per_word: 8,
            },
        },
        SimultaneousEventPolicyV1::NoSimultaneousEvents {
            witness: HybridWitnessIdV1::from_bytes(bytes(83)),
        },
        HybridTimeDomainV1 {
            start: 0.0,
            end: HybridTimeEndV1::Finite(4.0),
            event_cap: None,
        },
        CompactnessSemanticsV1::LocallyCompact {
            witness: HybridWitnessIdV1::from_bytes(bytes(84)),
        },
        AccumulationCandidateV1::None {
            no_claim: HybridNoClaimIdV1::from_bytes(bytes(85)),
        },
        ContinuationCategoryV1::Unique {
            rule: ContinuationRuleIdV1::from_bytes(bytes(86)),
            witness: HybridWitnessIdV1::from_bytes(bytes(87)),
        },
        HybridAnalysisBudgetV1 {
            max_event_word_len: 32,
            max_transitions: 1_000,
            max_wall_seconds: 1.0,
        },
    )
}

fn assert_problem_issue(ir: ZenoProblemIrV1, expected: HybridSemanticIssueV1) {
    with_cx(false, |cx| {
        let report = validate_zeno_problem_v1(ir, cx).expect_err("problem must refuse");
        assert!(
            report.issues().contains(&expected),
            "missing {expected:?} in {:?}",
            report.issues()
        );
    });
}

fn assert_claim_issue(
    draft: ZenoClaimDraftV1,
    problem: &ValidatedZenoProblemV1,
    regularized: Option<&ValidatedZenoProblemV1>,
    expected: HybridSemanticIssueV1,
) {
    with_cx(false, |cx| {
        let report = validate_zeno_claim_descriptor_v1(draft, problem, regularized, cx)
            .expect_err("claim descriptor must refuse");
        assert!(
            report.issues().contains(&expected),
            "missing {expected:?} in {:?}",
            report.issues()
        );
    });
}

fn set_valued_post() -> PostZenoStateV1 {
    PostZenoStateV1::SetValued {
        states: HybridStateSetIdV1::from_bytes(bytes(90)),
        rule: ContinuationRuleIdV1::from_bytes(bytes(52)),
        witness: HybridWitnessIdV1::from_bytes(bytes(91)),
    }
}

#[test]
fn relay_chatter_problem_replays_independent_of_mode_and_event_input_order() {
    let first = with_cx(false, |cx| {
        validate_zeno_problem_v1(relay_ir(), cx).unwrap()
    });
    assert!(first.has_zero_time_cycle());
    assert!(!first.has_nonunique_local_semantics());

    let mut reordered = relay_ir();
    reordered.modes.reverse();
    reordered.events.reverse();
    let reordered = with_cx(false, |cx| validate_zeno_problem_v1(reordered, cx).unwrap());
    assert_eq!(first.problem_id(), reordered.problem_id());
    assert_eq!(first.identity_receipt(), reordered.identity_receipt());

    let replay = with_cx(false, |cx| {
        validate_zeno_problem_v1(first.ir().clone(), cx).unwrap()
    });
    assert_eq!(first.identity_receipt(), replay.identity_receipt());

    let mut changed_version = relay_ir();
    changed_version.model_version = HybridModelVersionIdV1::from_bytes(bytes(92));
    let changed_version = with_cx(false, |cx| {
        validate_zeno_problem_v1(changed_version, cx).unwrap()
    });
    assert_ne!(first.problem_id(), changed_version.problem_id());
}

#[test]
fn event_cap_and_dense_trace_cannot_mint_a_zeno_theorem_descriptor() {
    let problem = with_cx(false, |cx| {
        validate_zeno_problem_v1(relay_ir(), cx).unwrap()
    });
    let capped = ZenoClaimDraftV1::new(
        problem.problem_id(),
        ZenoClassificationV1::CertifiedZeno {
            interval: HybridTimeIntervalV1 {
                earliest: 0.95,
                latest: 1.05,
            },
            states: HybridStateSetIdV1::from_bytes(bytes(90)),
            evidence: ZenoEvidenceReferenceV1::EventCap {
                trace: HybridEventTraceIdV1::from_bytes(bytes(51)),
            },
        },
        set_valued_post(),
    );
    assert_claim_issue(
        capped,
        &problem,
        None,
        HybridSemanticIssueV1::InsufficientTheoremEvidence,
    );

    let numerical = ZenoClaimDraftV1::new(
        problem.problem_id(),
        ZenoClassificationV1::CertifiedZeno {
            interval: HybridTimeIntervalV1 {
                earliest: 0.95,
                latest: 1.05,
            },
            states: HybridStateSetIdV1::from_bytes(bytes(90)),
            evidence: ZenoEvidenceReferenceV1::NumericalOnly {
                trace: HybridEventTraceIdV1::from_bytes(bytes(51)),
            },
        },
        set_valued_post(),
    );
    assert_claim_issue(
        numerical,
        &problem,
        None,
        HybridSemanticIssueV1::InsufficientTheoremEvidence,
    );

    let theorem_shaped = ZenoClaimDraftV1::new(
        problem.problem_id(),
        ZenoClassificationV1::CertifiedZeno {
            interval: HybridTimeIntervalV1 {
                earliest: 0.95,
                latest: 1.05,
            },
            states: HybridStateSetIdV1::from_bytes(bytes(90)),
            evidence: ZenoEvidenceReferenceV1::IntervalValidated {
                witness: HybridWitnessIdV1::from_bytes(bytes(93)),
            },
        },
        set_valued_post(),
    );
    let descriptor = with_cx(false, |cx| {
        validate_zeno_claim_descriptor_v1(theorem_shaped, &problem, None, cx).unwrap()
    });
    assert_eq!(
        descriptor.scientific_authority(),
        ZenoScientificAuthorityV1::ScientificCorrectnessNotProven
    );
}

#[test]
fn finite_separation_and_accumulation_follow_the_zero_dwell_event_graph() {
    let separated = with_cx(false, |cx| {
        validate_zeno_problem_v1(separated_ir(), cx).unwrap()
    });
    assert!(!separated.has_zero_time_cycle());
    let finite = ZenoClaimDraftV1::new(
        separated.problem_id(),
        ZenoClassificationV1::FiniteEventSeparation {
            minimum_separation: 0.25,
            evidence: ZenoEvidenceReferenceV1::Analytic {
                witness: HybridWitnessIdV1::from_bytes(bytes(94)),
            },
        },
        PostZenoStateV1::NotApplicable {
            justification: HybridNoClaimIdV1::from_bytes(bytes(95)),
        },
    );
    with_cx(false, |cx| {
        validate_zeno_claim_descriptor_v1(finite, &separated, None, cx).unwrap();
    });

    let impossible_zeno = ZenoClaimDraftV1::new(
        separated.problem_id(),
        ZenoClassificationV1::CertifiedZeno {
            interval: HybridTimeIntervalV1 {
                earliest: 1.0,
                latest: 1.0,
            },
            states: HybridStateSetIdV1::from_bytes(bytes(96)),
            evidence: ZenoEvidenceReferenceV1::Analytic {
                witness: HybridWitnessIdV1::from_bytes(bytes(97)),
            },
        },
        PostZenoStateV1::Unique {
            state: HybridStateIdV1::from_bytes(bytes(98)),
            rule: ContinuationRuleIdV1::from_bytes(bytes(86)),
            witness: HybridWitnessIdV1::from_bytes(bytes(87)),
        },
    );
    assert_claim_issue(
        impossible_zeno,
        &separated,
        None,
        HybridSemanticIssueV1::ZenoAccumulationCycleRequired,
    );

    let chatter = with_cx(false, |cx| {
        validate_zeno_problem_v1(relay_ir(), cx).unwrap()
    });
    let impossible_separation = ZenoClaimDraftV1::new(
        chatter.problem_id(),
        ZenoClassificationV1::FiniteEventSeparation {
            minimum_separation: 0.01,
            evidence: ZenoEvidenceReferenceV1::Analytic {
                witness: HybridWitnessIdV1::from_bytes(bytes(99)),
            },
        },
        PostZenoStateV1::NotApplicable {
            justification: HybridNoClaimIdV1::from_bytes(bytes(100)),
        },
    );
    assert_claim_issue(
        impossible_separation,
        &chatter,
        None,
        HybridSemanticIssueV1::FiniteSeparationContradictsEventGraph,
    );
}

#[test]
fn retained_continuation_category_distinguishes_unique_set_terminal_and_unresolved() {
    let mut unique_ir = relay_ir();
    let unique_rule = ContinuationRuleIdV1::from_bytes(bytes(130));
    let unique_witness = HybridWitnessIdV1::from_bytes(bytes(131));
    unique_ir.continuation = ContinuationCategoryV1::Unique {
        rule: unique_rule,
        witness: unique_witness,
    };
    let unique_problem = with_cx(false, |cx| validate_zeno_problem_v1(unique_ir, cx).unwrap());
    let unique_claim = ZenoClaimDraftV1::new(
        unique_problem.problem_id(),
        ZenoClassificationV1::CertifiedZeno {
            interval: HybridTimeIntervalV1 {
                earliest: 0.95,
                latest: 1.05,
            },
            states: HybridStateSetIdV1::from_bytes(bytes(132)),
            evidence: ZenoEvidenceReferenceV1::Analytic {
                witness: HybridWitnessIdV1::from_bytes(bytes(133)),
            },
        },
        PostZenoStateV1::Unique {
            state: HybridStateIdV1::from_bytes(bytes(134)),
            rule: unique_rule,
            witness: unique_witness,
        },
    );
    with_cx(false, |cx| {
        validate_zeno_claim_descriptor_v1(unique_claim, &unique_problem, None, cx).unwrap();
    });

    let mut terminal_ir = relay_ir();
    let terminal_rule = ContinuationRuleIdV1::from_bytes(bytes(135));
    terminal_ir.continuation = ContinuationCategoryV1::Terminal {
        rule: terminal_rule,
    };
    let terminal_problem = with_cx(false, |cx| {
        validate_zeno_problem_v1(terminal_ir, cx).unwrap()
    });
    let terminal_claim = ZenoClaimDraftV1::new(
        terminal_problem.problem_id(),
        ZenoClassificationV1::CertifiedZeno {
            interval: HybridTimeIntervalV1 {
                earliest: 0.95,
                latest: 1.05,
            },
            states: HybridStateSetIdV1::from_bytes(bytes(136)),
            evidence: ZenoEvidenceReferenceV1::IntervalValidated {
                witness: HybridWitnessIdV1::from_bytes(bytes(137)),
            },
        },
        PostZenoStateV1::Terminal {
            rule: terminal_rule,
            witness: HybridWitnessIdV1::from_bytes(bytes(138)),
        },
    );
    with_cx(false, |cx| {
        validate_zeno_claim_descriptor_v1(terminal_claim, &terminal_problem, None, cx).unwrap();
    });

    assert_ne!(unique_problem.problem_id(), terminal_problem.problem_id());
}

fn simultaneous_ir() -> ZenoProblemIrV1 {
    let mut ir = separated_ir();
    let group = SimultaneityGroupIdV1::from_bytes(bytes(101));
    ir.events[0].simultaneity = EventSimultaneityV1::Group { group };
    let mut second = ir.events[0].clone();
    second.event = HybridEventIdV1::from_bytes(bytes(102));
    second.guard = HybridGuardIdV1::from_bytes(bytes(103));
    second.reset = ResetSemanticsV1::SetValued {
        relation: ResetRelationIdV1::from_bytes(bytes(104)),
        targets: vec![ir.modes[0].mode, ir.modes[1].mode],
        states: HybridStateSetIdV1::from_bytes(bytes(105)),
    };
    ir.events.push(second);
    ir.simultaneous_policy = SimultaneousEventPolicyV1::SetValued {
        outcomes: HybridStateSetIdV1::from_bytes(bytes(106)),
    };
    ir.continuation = ContinuationCategoryV1::SetValued {
        rule: ContinuationRuleIdV1::from_bytes(bytes(107)),
    };
    ir
}

#[test]
fn simultaneous_grazing_and_nonunique_resets_never_become_unique_silently() {
    let simultaneous = with_cx(false, |cx| {
        validate_zeno_problem_v1(simultaneous_ir(), cx).unwrap()
    });
    assert!(simultaneous.has_nonunique_local_semantics());

    let mut unique_overclaim = simultaneous_ir();
    unique_overclaim.continuation = ContinuationCategoryV1::Unique {
        rule: ContinuationRuleIdV1::from_bytes(bytes(108)),
        witness: HybridWitnessIdV1::from_bytes(bytes(109)),
    };
    assert_problem_issue(
        unique_overclaim,
        HybridSemanticIssueV1::UniqueContinuationUnsupported,
    );

    let mut grazing = separated_ir();
    grazing.events[0].crossing = CrossingSemanticsV1::Grazing {
        witness: HybridWitnessIdV1::from_bytes(bytes(110)),
    };
    assert_problem_issue(
        grazing,
        HybridSemanticIssueV1::UniqueContinuationUnsupported,
    );

    let wrong_post = ZenoClaimDraftV1::new(
        simultaneous.problem_id(),
        ZenoClassificationV1::CertifiedZeno {
            interval: HybridTimeIntervalV1 {
                earliest: 1.0,
                latest: 1.1,
            },
            states: HybridStateSetIdV1::from_bytes(bytes(111)),
            evidence: ZenoEvidenceReferenceV1::Analytic {
                witness: HybridWitnessIdV1::from_bytes(bytes(112)),
            },
        },
        PostZenoStateV1::Unique {
            state: HybridStateIdV1::from_bytes(bytes(113)),
            rule: ContinuationRuleIdV1::from_bytes(bytes(108)),
            witness: HybridWitnessIdV1::from_bytes(bytes(109)),
        },
    );
    assert_claim_issue(
        wrong_post,
        &simultaneous,
        None,
        HybridSemanticIssueV1::PostZenoSemanticsMismatch,
    );

    let mut invalid_priority = simultaneous_ir();
    invalid_priority.simultaneous_policy = SimultaneousEventPolicyV1::TotalPriority {
        ordered_events: vec![invalid_priority.events[0].event],
        witness: HybridWitnessIdV1::from_bytes(bytes(114)),
    };
    assert_problem_issue(
        invalid_priority,
        HybridSemanticIssueV1::InvalidPriorityOrder,
    );

    let mut singleton = separated_ir();
    singleton.events[0].simultaneity = EventSimultaneityV1::Group {
        group: SimultaneityGroupIdV1::from_bytes(bytes(115)),
    };
    singleton.simultaneous_policy = SimultaneousEventPolicyV1::SetValued {
        outcomes: HybridStateSetIdV1::from_bytes(bytes(116)),
    };
    singleton.continuation = ContinuationCategoryV1::SetValued {
        rule: ContinuationRuleIdV1::from_bytes(bytes(117)),
    };
    assert_problem_issue(singleton, HybridSemanticIssueV1::SingletonSimultaneityGroup);
}

#[test]
fn compliant_regularization_has_distinct_lineage_and_no_silent_equivalence() {
    let original_ir = relay_ir();
    let source_model = original_ir.model;
    let source_version = original_ir.model_version;
    let original = with_cx(false, |cx| {
        validate_zeno_problem_v1(original_ir, cx).unwrap()
    });
    let regularization = HybridRegularizationIdV1::from_bytes(bytes(118));
    let no_equivalence = HybridNoClaimIdV1::from_bytes(bytes(119));
    let mut regularized_ir = relay_ir();
    regularized_ir.model = HybridModelIdV1::from_bytes(bytes(120));
    regularized_ir.model_version = HybridModelVersionIdV1::from_bytes(bytes(121));
    regularized_ir.lineage = HybridModelLineageV1::Regularized {
        source_model,
        source_version,
        regularization,
        no_equivalence,
    };
    let regularized = with_cx(false, |cx| {
        validate_zeno_problem_v1(regularized_ir, cx).unwrap()
    });
    assert_ne!(original.problem_id(), regularized.problem_id());

    let claim = ZenoClaimDraftV1::new(
        original.problem_id(),
        ZenoClassificationV1::RegularizedModel {
            regularized_problem: regularized.problem_id(),
            regularization,
            no_equivalence,
        },
        PostZenoStateV1::Unresolved {
            no_claim: HybridNoClaimIdV1::from_bytes(bytes(122)),
        },
    );
    let admitted = with_cx(false, |cx| {
        validate_zeno_claim_descriptor_v1(claim, &original, Some(&regularized), cx).unwrap()
    });
    assert_eq!(
        admitted.scientific_authority(),
        ZenoScientificAuthorityV1::ScientificCorrectnessNotProven
    );
    assert_claim_issue(
        claim,
        &original,
        None,
        HybridSemanticIssueV1::RegularizationLineageMismatch,
    );

    let mut self_reference = relay_ir();
    self_reference.lineage = HybridModelLineageV1::Regularized {
        source_model: self_reference.model,
        source_version: self_reference.model_version,
        regularization,
        no_equivalence,
    };
    assert_problem_issue(
        self_reference,
        HybridSemanticIssueV1::RegularizationSelfReference,
    );
}

#[test]
fn numerical_warning_stays_unknown_after_the_event_window() {
    let problem = with_cx(false, |cx| {
        validate_zeno_problem_v1(relay_ir(), cx).unwrap()
    });
    let warning = ZenoClaimDraftV1::new(
        problem.problem_id(),
        ZenoClassificationV1::NumericalEventDensityWarning {
            trace: HybridEventTraceIdV1::from_bytes(bytes(123)),
            observed_events: 1000,
            window: 0.001,
        },
        PostZenoStateV1::Unresolved {
            no_claim: HybridNoClaimIdV1::from_bytes(bytes(124)),
        },
    );
    let first = with_cx(false, |cx| {
        validate_zeno_claim_descriptor_v1(warning, &problem, None, cx).unwrap()
    });
    let replay = with_cx(false, |cx| {
        validate_zeno_claim_descriptor_v1(*first.draft(), &problem, None, cx).unwrap()
    });
    assert_eq!(first.identity_receipt(), replay.identity_receipt());

    let overclaim = ZenoClaimDraftV1::new(
        problem.problem_id(),
        warning.classification(),
        set_valued_post(),
    );
    assert_claim_issue(
        overclaim,
        &problem,
        None,
        HybridSemanticIssueV1::PostZenoSemanticsMismatch,
    );
}

fn with_problem_schema(ir: ZenoProblemIrV1, schema_version: u32) -> ZenoProblemIrV1 {
    ZenoProblemIrV1::with_schema_version(
        schema_version,
        ir.model,
        ir.model_version,
        ir.lineage,
        ir.frame,
        ir.state_units,
        ir.initial_states,
        ir.time_scale,
        ir.modes,
        ir.events,
        ir.event_language,
        ir.simultaneous_policy,
        ir.time_domain,
        ir.compactness,
        ir.accumulation_candidate,
        ir.continuation,
        ir.budget,
    )
}

#[test]
fn invalid_references_degrees_windows_versions_and_caps_refuse() {
    let mut missing_mode = separated_ir();
    missing_mode.events[0].source_mode = HybridModeIdV1::from_bytes(bytes(200));
    assert_problem_issue(missing_mode, HybridSemanticIssueV1::UnknownModeReference);

    let mut bad_dae = separated_ir();
    bad_dae.modes[0].class = ContinuousDynamicsClassV1::AdmittedDae {
        index: 0,
        constraint: DaeConstraintIdV1::from_bytes(bytes(201)),
    };
    assert_problem_issue(
        bad_dae,
        HybridSemanticIssueV1::InvalidValue {
            field: HybridFieldV1::DaeIndex,
        },
    );

    let mut outside = separated_ir();
    outside.accumulation_candidate = AccumulationCandidateV1::Window {
        earliest: 3.0,
        latest: 5.0,
        states: HybridStateSetIdV1::from_bytes(bytes(202)),
        trace: HybridEventTraceIdV1::from_bytes(bytes(203)),
    };
    assert_problem_issue(
        outside,
        HybridSemanticIssueV1::AccumulationOutsideTimeDomain,
    );

    assert_problem_issue(
        with_problem_schema(separated_ir(), 99),
        HybridSemanticIssueV1::UnsupportedSchemaVersion {
            found: 99,
            supported: ZENO_PROBLEM_SCHEMA_VERSION_V1,
        },
    );

    let template = mode(210);
    let mut oversized = separated_ir();
    oversized.modes = vec![template; MAX_HYBRID_MODES_V1 + 1];
    assert_problem_issue(
        oversized,
        HybridSemanticIssueV1::TooMany {
            collection: HybridCollectionV1::Modes,
            found: MAX_HYBRID_MODES_V1 + 1,
            limit: MAX_HYBRID_MODES_V1,
        },
    );

    let mut oversized_events = separated_ir();
    let target = oversized_events.modes[0].mode;
    let mut oversized_reset = oversized_events.events[0].clone();
    oversized_reset.reset = ResetSemanticsV1::SetValued {
        relation: ResetRelationIdV1::from_bytes(bytes(204)),
        targets: vec![target; MAX_RESET_TARGETS_V1 + 1],
        states: HybridStateSetIdV1::from_bytes(bytes(205)),
    };
    oversized_events.events = vec![oversized_events.events[0].clone(); MAX_HYBRID_EVENTS_V1 + 1];
    oversized_events.events[0] = oversized_reset;
    let report = with_cx(false, |cx| {
        validate_zeno_problem_v1(oversized_events, cx).expect_err("top-level cap refuses")
    });
    assert_eq!(
        report.issues(),
        &[HybridSemanticIssueV1::TooMany {
            collection: HybridCollectionV1::Events,
            found: MAX_HYBRID_EVENTS_V1 + 1,
            limit: MAX_HYBRID_EVENTS_V1,
        }]
    );

    let mut oversized_priority = separated_ir();
    oversized_priority.simultaneous_policy = SimultaneousEventPolicyV1::TotalPriority {
        ordered_events: vec![oversized_priority.events[0].event; MAX_HYBRID_EVENTS_V1 + 1],
        witness: HybridWitnessIdV1::from_bytes(bytes(206)),
    };
    let report = with_cx(false, |cx| {
        validate_zeno_problem_v1(oversized_priority, cx).expect_err("priority cap refuses")
    });
    assert_eq!(
        report.issues(),
        &[HybridSemanticIssueV1::TooMany {
            collection: HybridCollectionV1::PriorityEvents,
            found: MAX_HYBRID_EVENTS_V1 + 1,
            limit: MAX_HYBRID_EVENTS_V1,
        }]
    );
}

#[test]
fn reset_target_cap_report_is_permutation_stable() {
    let mut first = separated_ir();
    let source = first.modes[0].mode;
    let target = first.modes[1].mode;
    let mut smaller = event(
        230,
        source,
        target,
        InteractionLawV1::None {
            justification: HybridNoClaimIdV1::from_bytes(bytes(231)),
        },
        DwellSemanticsV1::ZeroAllowed,
    );
    smaller.reset = ResetSemanticsV1::SetValued {
        relation: ResetRelationIdV1::from_bytes(bytes(232)),
        targets: vec![target; MAX_RESET_TARGETS_V1 + 1],
        states: HybridStateSetIdV1::from_bytes(bytes(233)),
    };
    let mut larger = event(
        240,
        source,
        target,
        InteractionLawV1::None {
            justification: HybridNoClaimIdV1::from_bytes(bytes(241)),
        },
        DwellSemanticsV1::ZeroAllowed,
    );
    larger.reset = ResetSemanticsV1::SetValued {
        relation: ResetRelationIdV1::from_bytes(bytes(242)),
        targets: vec![target; MAX_RESET_TARGETS_V1 + 2],
        states: HybridStateSetIdV1::from_bytes(bytes(243)),
    };
    first.events = vec![smaller, larger];
    let mut reversed = first.clone();
    reversed.events.reverse();

    let first_report = with_cx(false, |cx| {
        validate_zeno_problem_v1(first, cx).expect_err("reset cap refuses")
    });
    let reversed_report = with_cx(false, |cx| {
        validate_zeno_problem_v1(reversed, cx).expect_err("reset cap refuses")
    });
    assert_eq!(first_report, reversed_report);
    assert_eq!(
        first_report.issues(),
        &[HybridSemanticIssueV1::TooMany {
            collection: HybridCollectionV1::ResetTargets,
            found: MAX_RESET_TARGETS_V1 + 2,
            limit: MAX_RESET_TARGETS_V1,
        }]
    );
}

#[test]
fn total_reset_target_cap_report_is_permutation_stable() {
    let mut first = separated_ir();
    let target = first.modes[1].mode;
    let template = first.events[0].clone();
    let event_count = MAX_TOTAL_RESET_TARGETS_V1 / MAX_RESET_TARGETS_V1 + 1;
    first.events = (0..event_count)
        .map(|index| {
            let seed = u8::try_from(index + 1).expect("fixture event count fits u8");
            let mut event = template.clone();
            event.event = HybridEventIdV1::from_bytes(bytes(seed));
            event.reset = ResetSemanticsV1::SetValued {
                relation: ResetRelationIdV1::from_bytes(bytes(180 + seed)),
                targets: vec![target; MAX_RESET_TARGETS_V1],
                states: HybridStateSetIdV1::from_bytes(bytes(210 + seed)),
            };
            event
        })
        .collect();
    let mut reversed = first.clone();
    reversed.events.reverse();

    let first_report = with_cx(false, |cx| {
        validate_zeno_problem_v1(first, cx).expect_err("total reset cap refuses")
    });
    let reversed_report = with_cx(false, |cx| {
        validate_zeno_problem_v1(reversed, cx).expect_err("total reset cap refuses")
    });
    assert_eq!(first_report, reversed_report);
    assert_eq!(
        first_report.issues(),
        &[HybridSemanticIssueV1::TooMany {
            collection: HybridCollectionV1::TotalResetTargets,
            found: event_count * MAX_RESET_TARGETS_V1,
            limit: MAX_TOTAL_RESET_TARGETS_V1,
        }]
    );
}

#[test]
fn pre_cancelled_validation_publishes_no_problem_or_claim_identity() {
    with_cx(true, |cx| {
        let report = validate_zeno_problem_v1(separated_ir(), cx).unwrap_err();
        assert_eq!(report.issues(), &[HybridSemanticIssueV1::Cancelled]);
    });

    let problem = with_cx(false, |cx| {
        validate_zeno_problem_v1(separated_ir(), cx).unwrap()
    });
    let draft = ZenoClaimDraftV1::new(
        problem.problem_id(),
        ZenoClassificationV1::Unknown {
            no_claim: HybridNoClaimIdV1::from_bytes(bytes(220)),
        },
        PostZenoStateV1::Unresolved {
            no_claim: HybridNoClaimIdV1::from_bytes(bytes(221)),
        },
    );
    with_cx(true, |cx| {
        let report = validate_zeno_claim_descriptor_v1(draft, &problem, None, cx).unwrap_err();
        assert_eq!(report.issues(), &[HybridSemanticIssueV1::Cancelled]);
    });
}
