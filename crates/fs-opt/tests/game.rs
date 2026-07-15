//! G0/G3/G4/G5 conformance for RE.Q1 quantified-game semantics.

#![allow(clippy::wildcard_imports)]
#![allow(
    clippy::too_many_lines,
    reason = "each long RE.Q1 case keeps one information or quantifier refusal narrative explicit"
)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_opt::game::*;

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
                seed: 0x4741_4D45,
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

fn context() -> StateSetContextV1 {
    StateSetContextV1 {
        state_space: GameStateSpaceIdV1::from_bytes(bytes(3)),
        frame: GameFrameIdV1::from_bytes(bytes(4)),
        units: GameUnitSystemIdV1::from_bytes(bytes(5)),
        model_version: GameModelVersionIdV1::from_bytes(bytes(2)),
    }
}

fn base_ir() -> GameProblemIrV1 {
    let model = GameModelSpecV1 {
        model: GameModelIdV1::from_bytes(bytes(1)),
        version: GameModelVersionIdV1::from_bytes(bytes(2)),
        state_space: GameStateSpaceIdV1::from_bytes(bytes(3)),
        frame: GameFrameIdV1::from_bytes(bytes(4)),
        state_units: GameUnitSystemIdV1::from_bytes(bytes(5)),
        class: GameModelClassV1::DifferentialGame,
    };
    let initial = InitialSetV1 {
        set: InitialSetIdV1::from_bytes(bytes(6)),
        context: context(),
    };
    let target = TargetSetV1 {
        set: TargetSetIdV1::from_bytes(bytes(7)),
        context: context(),
    };
    let unsafe_set = UnsafeSetV1 {
        set: UnsafeSetIdV1::from_bytes(bytes(8)),
        context: context(),
    };
    let controls = ControlSetV1 {
        set: ControlSetIdV1::from_bytes(bytes(9)),
        units: GameUnitSystemIdV1::from_bytes(bytes(10)),
        model_version: model.version,
    };
    let disturbances = DisturbanceSetV1 {
        set: DisturbanceSetIdV1::from_bytes(bytes(11)),
        units: GameUnitSystemIdV1::from_bytes(bytes(12)),
        model_version: model.version,
    };
    let parameters = ParameterSetV1 {
        set: ParameterSetIdV1::from_bytes(bytes(13)),
        units: GameUnitSystemIdV1::from_bytes(bytes(14)),
        model_version: model.version,
    };
    let quantifiers = GameQuantifierPrefixV1::new(vec![
        GameQuantifierClauseV1 {
            quantifier: GameQuantifierV1::Exists,
            variable: GameVariableV1::Control,
            domain: GameQuantifierDomainV1::Control(controls.set),
        },
        GameQuantifierClauseV1 {
            quantifier: GameQuantifierV1::ForAll,
            variable: GameVariableV1::Disturbance,
            domain: GameQuantifierDomainV1::Disturbance(disturbances.set),
        },
        GameQuantifierClauseV1 {
            quantifier: GameQuantifierV1::ForAll,
            variable: GameVariableV1::Parameter,
            domain: GameQuantifierDomainV1::Parameter(parameters.set),
        },
    ]);
    let information = InformationPatternV1 {
        grants: vec![
            InformationGrantV1 {
                player: GamePlayerV1::Controller,
                subject: InformationSubjectV1::Parameter,
                observation: GameObservationMapIdV1::from_bytes(bytes(15)),
                availability: ObservationAvailabilityV1::InitialOnly,
            },
            InformationGrantV1 {
                player: GamePlayerV1::Controller,
                subject: InformationSubjectV1::Disturbance,
                observation: GameObservationMapIdV1::from_bytes(bytes(16)),
                availability: ObservationAvailabilityV1::Hidden,
            },
            InformationGrantV1 {
                player: GamePlayerV1::Controller,
                subject: InformationSubjectV1::State,
                observation: GameObservationMapIdV1::from_bytes(bytes(17)),
                availability: ObservationAvailabilityV1::Current,
            },
        ],
    };
    let strategies = vec![GameStrategySpecV1 {
        player: GamePlayerV1::Controller,
        representation: StrategyRepresentationV1::StateFeedback {
            artifact: GameStrategyArtifactIdV1::from_bytes(bytes(18)),
        },
        dependencies: vec![
            StrategyDependencyV1 {
                subject: InformationSubjectV1::State,
                access: StrategyAccessV1::Current,
            },
            StrategyDependencyV1 {
                subject: InformationSubjectV1::Parameter,
                access: StrategyAccessV1::InitialOnly,
            },
        ],
    }];
    GameProblemIrV1::new(
        model,
        initial,
        target,
        unsafe_set,
        controls,
        disturbances,
        parameters,
        GameClaimSpecV1 {
            objective: GameObjectiveV1::ReachAvoid,
            polarity: GameProofPolarityV1::Inner,
        },
        quantifiers,
        information,
        strategies,
        GameHorizonV1::Finite {
            start: 0.0,
            end: 10.0,
            unit: GameTimeUnitIdV1::from_bytes(bytes(19)),
            seconds_per_unit: 1.0,
        },
        GameStoppingSemanticsV1::FirstTargetOrUnsafe,
        GameCompositionV1::Atomic,
        GameAnalysisBudgetV1 {
            cell_decomposition: GameCellDecompositionIdV1::from_bytes(bytes(47)),
            max_cells: 100_000,
            max_transitions: 1_000_000,
            max_strategy_nodes: 10_000,
            max_wall_seconds: 5.0,
        },
    )
}

fn validate(ir: GameProblemIrV1) -> ValidatedGameProblemV1 {
    with_cx(false, |cx| {
        validate_game_problem_v1(ir, cx).expect("fixture admits")
    })
}

fn assert_issue(ir: GameProblemIrV1, expected: GameSemanticIssueV1) {
    let report = with_cx(false, |cx| {
        validate_game_problem_v1(ir, cx).expect_err("fixture must refuse")
    });
    assert!(
        report.issues().contains(&expected),
        "missing {expected:?} in {:?}",
        report.issues()
    );
}

fn controller_strategy(ir: &mut GameProblemIrV1) -> &mut GameStrategySpecV1 {
    ir.strategies
        .iter_mut()
        .find(|strategy| strategy.player == GamePlayerV1::Controller)
        .expect("controller fixture")
}

fn set_disturbance_grant(ir: &mut GameProblemIrV1, availability: ObservationAvailabilityV1) {
    ir.information
        .grants
        .iter_mut()
        .find(|grant| {
            grant.player == GamePlayerV1::Controller
                && grant.subject == InformationSubjectV1::Disturbance
        })
        .expect("disturbance grant")
        .availability = availability;
}

fn swap_prefix(ir: &mut GameProblemIrV1) {
    let clauses = [
        ir.quantifiers.clauses()[0],
        ir.quantifiers.clauses()[1],
        ir.quantifiers.clauses()[2],
    ];
    ir.quantifiers = GameQuantifierPrefixV1::new(vec![clauses[1], clauses[0], clauses[2]]);
}

#[test]
fn canonical_replay_is_order_free_except_for_quantifiers() {
    let first = validate(base_ir());
    let mut reordered = base_ir();
    reordered.information.grants.reverse();
    reordered.strategies[0].dependencies.reverse();
    let second = validate(reordered);
    assert_eq!(first.identity_receipt(), second.identity_receipt());
    assert_eq!(
        first.claim_availability(),
        GameClaimAvailabilityV1::Eligible {
            polarity: GameProofPolarityV1::Inner,
        }
    );

    let mut swapped = base_ir();
    swap_prefix(&mut swapped);
    let rendered = swapped.quantifiers.to_string();
    assert!(rendered.starts_with("forall disturbance in disturbance-set("));
    assert!(rendered.contains("exists control in control-set("));
    assert!(rendered.contains("disturbance-set(0b0b0b"));
    assert!(rendered.contains("control-set(090909"));
    let summary = swapped.to_string();
    assert!(summary.contains("quantifiers=[forall disturbance"));
    assert!(summary.contains("controller:state-feedback"));
    let swapped = validate(swapped);
    assert_ne!(first.problem_id(), swapped.problem_id());
}

#[test]
fn open_loop_feedback_and_prefix_order_are_not_confusable() {
    let feedback = validate(base_ir());
    let mut open_loop = base_ir();
    let strategy = controller_strategy(&mut open_loop);
    strategy.representation = StrategyRepresentationV1::OpenLoop {
        artifact: GameStrategyArtifactIdV1::from_bytes(bytes(20)),
    };
    strategy
        .dependencies
        .retain(|dependency| dependency.subject == InformationSubjectV1::Parameter);
    let open_loop = validate(open_loop);
    assert_ne!(feedback.problem_id(), open_loop.problem_id());

    let mut clairvoyant = base_ir();
    swap_prefix(&mut clairvoyant);
    let strategy = controller_strategy(&mut clairvoyant);
    strategy.representation = StrategyRepresentationV1::OpenLoop {
        artifact: GameStrategyArtifactIdV1::from_bytes(bytes(21)),
    };
    strategy.dependencies.clear();
    assert_issue(
        clairvoyant,
        GameSemanticIssueV1::ClairvoyantQuantifierLowering,
    );

    let mut causal = base_ir();
    swap_prefix(&mut causal);
    set_disturbance_grant(&mut causal, ObservationAvailabilityV1::Current);
    let strategy = controller_strategy(&mut causal);
    strategy.representation = StrategyRepresentationV1::NonanticipativeFeedback {
        artifact: GameStrategyArtifactIdV1::from_bytes(bytes(22)),
        memory_states: 32,
    };
    strategy.dependencies.push(StrategyDependencyV1 {
        subject: InformationSubjectV1::Disturbance,
        access: StrategyAccessV1::Current,
    });
    validate(causal);
}

#[test]
fn dependency_free_open_loop_needs_no_dummy_information_grant() {
    let mut open_loop = base_ir();
    open_loop.information.grants.clear();
    let strategy = controller_strategy(&mut open_loop);
    strategy.representation = StrategyRepresentationV1::OpenLoop {
        artifact: GameStrategyArtifactIdV1::from_bytes(bytes(50)),
    };
    strategy.dependencies.clear();

    validate(open_loop);
}

#[test]
fn nonanticipative_memory_must_fit_the_strategy_node_budget() {
    let mut excessive = base_ir();
    let limit = excessive.budget.max_strategy_nodes;
    controller_strategy(&mut excessive).representation =
        StrategyRepresentationV1::NonanticipativeFeedback {
            artifact: GameStrategyArtifactIdV1::from_bytes(bytes(51)),
            memory_states: u32::try_from(limit + 1).expect("fixture limit fits u32"),
        };
    assert_issue(
        excessive,
        GameSemanticIssueV1::InvalidValue {
            field: GameFieldV1::StrategyMemory,
        },
    );

    let mut boundary = base_ir();
    controller_strategy(&mut boundary).representation =
        StrategyRepresentationV1::NonanticipativeFeedback {
            artifact: GameStrategyArtifactIdV1::from_bytes(bytes(52)),
            memory_states: u32::try_from(limit).expect("fixture limit fits u32"),
        };
    validate(boundary);
}

#[test]
fn hidden_and_delayed_disturbances_enforce_information_timing() {
    let mut hidden = base_ir();
    let strategy = controller_strategy(&mut hidden);
    strategy.representation = StrategyRepresentationV1::NonanticipativeFeedback {
        artifact: GameStrategyArtifactIdV1::from_bytes(bytes(23)),
        memory_states: 8,
    };
    strategy.dependencies.push(StrategyDependencyV1 {
        subject: InformationSubjectV1::Disturbance,
        access: StrategyAccessV1::Current,
    });
    assert_issue(hidden, GameSemanticIssueV1::HiddenOrMissingObservation);

    let mut too_early = base_ir();
    set_disturbance_grant(
        &mut too_early,
        ObservationAvailabilityV1::Delayed { lag: 0.5 },
    );
    controller_strategy(&mut too_early)
        .dependencies
        .push(StrategyDependencyV1 {
            subject: InformationSubjectV1::Disturbance,
            access: StrategyAccessV1::Delayed { lag: 0.25 },
        });
    assert_issue(too_early, GameSemanticIssueV1::ObservationTimingMismatch);

    let mut unavailable_at_initial_choice = base_ir();
    set_disturbance_grant(
        &mut unavailable_at_initial_choice,
        ObservationAvailabilityV1::Delayed { lag: 0.5 },
    );
    let strategy = controller_strategy(&mut unavailable_at_initial_choice);
    strategy.representation = StrategyRepresentationV1::OpenLoop {
        artifact: GameStrategyArtifactIdV1::from_bytes(bytes(48)),
    };
    strategy.dependencies = vec![StrategyDependencyV1 {
        subject: InformationSubjectV1::Disturbance,
        access: StrategyAccessV1::InitialOnly,
    }];
    assert_issue(
        unavailable_at_initial_choice,
        GameSemanticIssueV1::ObservationTimingMismatch,
    );

    let mut positive_lag_current = base_ir();
    set_disturbance_grant(
        &mut positive_lag_current,
        ObservationAvailabilityV1::Delayed { lag: 0.5 },
    );
    let strategy = controller_strategy(&mut positive_lag_current);
    strategy.representation = StrategyRepresentationV1::NonanticipativeFeedback {
        artifact: GameStrategyArtifactIdV1::from_bytes(bytes(56)),
        memory_states: 8,
    };
    strategy.dependencies.push(StrategyDependencyV1 {
        subject: InformationSubjectV1::Disturbance,
        access: StrategyAccessV1::Current,
    });
    assert_issue(
        positive_lag_current,
        GameSemanticIssueV1::ObservationTimingMismatch,
    );

    let mut zero_lag_initial = base_ir();
    set_disturbance_grant(
        &mut zero_lag_initial,
        ObservationAvailabilityV1::Delayed { lag: 0.0 },
    );
    let strategy = controller_strategy(&mut zero_lag_initial);
    strategy.representation = StrategyRepresentationV1::OpenLoop {
        artifact: GameStrategyArtifactIdV1::from_bytes(bytes(49)),
    };
    strategy.dependencies = vec![StrategyDependencyV1 {
        subject: InformationSubjectV1::Disturbance,
        access: StrategyAccessV1::InitialOnly,
    }];
    validate(zero_lag_initial);

    let mut zero_lag_current = base_ir();
    set_disturbance_grant(
        &mut zero_lag_current,
        ObservationAvailabilityV1::Delayed { lag: 0.0 },
    );
    let strategy = controller_strategy(&mut zero_lag_current);
    strategy.representation = StrategyRepresentationV1::NonanticipativeFeedback {
        artifact: GameStrategyArtifactIdV1::from_bytes(bytes(53)),
        memory_states: 8,
    };
    strategy.dependencies.push(StrategyDependencyV1 {
        subject: InformationSubjectV1::Disturbance,
        access: StrategyAccessV1::Current,
    });
    validate(zero_lag_current);

    let mut zero_lag_history = base_ir();
    set_disturbance_grant(
        &mut zero_lag_history,
        ObservationAvailabilityV1::Delayed { lag: 0.0 },
    );
    let strategy = controller_strategy(&mut zero_lag_history);
    strategy.representation = StrategyRepresentationV1::NonanticipativeFeedback {
        artifact: GameStrategyArtifactIdV1::from_bytes(bytes(57)),
        memory_states: 8,
    };
    strategy.dependencies.push(StrategyDependencyV1 {
        subject: InformationSubjectV1::Disturbance,
        access: StrategyAccessV1::HistoryThroughCurrent,
    });
    assert_issue(
        zero_lag_history,
        GameSemanticIssueV1::ObservationTimingMismatch,
    );

    let mut admitted = base_ir();
    set_disturbance_grant(
        &mut admitted,
        ObservationAvailabilityV1::Delayed { lag: 0.5 },
    );
    let strategy = controller_strategy(&mut admitted);
    strategy.representation = StrategyRepresentationV1::NonanticipativeFeedback {
        artifact: GameStrategyArtifactIdV1::from_bytes(bytes(24)),
        memory_states: 8,
    };
    strategy.dependencies.push(StrategyDependencyV1 {
        subject: InformationSubjectV1::Disturbance,
        access: StrategyAccessV1::Delayed { lag: 0.5 },
    });
    validate(admitted);
}

#[test]
fn controller_and_disturbance_quantifiers_match_player_roles() {
    let mut universal_control = base_ir();
    let clauses = [
        universal_control.quantifiers.clauses()[0],
        universal_control.quantifiers.clauses()[1],
        universal_control.quantifiers.clauses()[2],
    ];
    universal_control.quantifiers = GameQuantifierPrefixV1::new(vec![
        GameQuantifierClauseV1 {
            quantifier: GameQuantifierV1::ForAll,
            ..clauses[0]
        },
        clauses[1],
        clauses[2],
    ]);
    assert_issue(
        universal_control,
        GameSemanticIssueV1::QuantifierPolarityMismatch {
            variable: GameVariableV1::Control,
            found: GameQuantifierV1::ForAll,
            required: GameQuantifierV1::Exists,
        },
    );

    let mut existential_disturbance = base_ir();
    let clauses = [
        existential_disturbance.quantifiers.clauses()[0],
        existential_disturbance.quantifiers.clauses()[1],
        existential_disturbance.quantifiers.clauses()[2],
    ];
    existential_disturbance.quantifiers = GameQuantifierPrefixV1::new(vec![
        clauses[0],
        GameQuantifierClauseV1 {
            quantifier: GameQuantifierV1::Exists,
            ..clauses[1]
        },
        clauses[2],
    ]);
    assert_issue(
        existential_disturbance,
        GameSemanticIssueV1::QuantifierPolarityMismatch {
            variable: GameVariableV1::Disturbance,
            found: GameQuantifierV1::Exists,
            required: GameQuantifierV1::ForAll,
        },
    );
}

#[test]
fn repeated_quantifier_diagnostics_are_canonical() {
    let mut repeated = base_ir();
    let clauses = [
        repeated.quantifiers.clauses()[0],
        repeated.quantifiers.clauses()[1],
        repeated.quantifiers.clauses()[2],
    ];
    repeated.quantifiers = GameQuantifierPrefixV1::new(vec![
        clauses[0], clauses[0], clauses[0], clauses[1], clauses[2],
    ]);

    let report = with_cx(false, |cx| {
        validate_game_problem_v1(repeated, cx).expect_err("duplicate control quantifiers refuse")
    });
    assert_eq!(
        report.issues(),
        &[GameSemanticIssueV1::DuplicateQuantifiedVariable {
            variable: GameVariableV1::Control,
        }]
    );
}

#[test]
fn infinite_horizon_is_admitted_only_as_unknown() {
    let mut infinite = base_ir();
    infinite.horizon = GameHorizonV1::Infinite {
        start: 0.0,
        unit: GameTimeUnitIdV1::from_bytes(bytes(25)),
        seconds_per_unit: 1.0,
        no_claim: GameNoClaimIdV1::from_bytes(bytes(26)),
    };
    infinite.stopping = GameStoppingSemanticsV1::Never;
    let admitted = validate(infinite);
    assert_eq!(
        admitted.claim_availability(),
        GameClaimAvailabilityV1::Unknown {
            reason: GameUnknownReasonV1::InfiniteHorizon,
        }
    );

    let mut invalid = base_ir();
    invalid.horizon = GameHorizonV1::Infinite {
        start: 0.0,
        unit: GameTimeUnitIdV1::from_bytes(bytes(25)),
        seconds_per_unit: 1.0,
        no_claim: GameNoClaimIdV1::from_bytes(bytes(26)),
    };
    invalid.stopping = GameStoppingSemanticsV1::FixedHorizon;
    assert_issue(invalid, GameSemanticIssueV1::StoppingSemanticsMismatch);
}

#[test]
fn objective_and_inner_outer_unknown_polarities_remain_distinct() {
    let inner = validate(base_ir());
    let mut outer = base_ir();
    outer.claim.objective = GameObjectiveV1::Viability;
    outer.claim.polarity = GameProofPolarityV1::Outer;
    outer.stopping = GameStoppingSemanticsV1::FixedHorizon;
    let outer = validate(outer);
    assert_ne!(inner.problem_id(), outer.problem_id());
    assert_eq!(
        outer.claim_availability(),
        GameClaimAvailabilityV1::Eligible {
            polarity: GameProofPolarityV1::Outer,
        }
    );

    let mut unknown = base_ir();
    unknown.claim.polarity = GameProofPolarityV1::Unknown {
        no_claim: GameNoClaimIdV1::from_bytes(bytes(46)),
    };
    assert_eq!(
        validate(unknown).claim_availability(),
        GameClaimAvailabilityV1::Unknown {
            reason: GameUnknownReasonV1::PolarityUnresolved,
        }
    );
}

#[test]
fn stopping_semantics_preserve_the_full_objective_obligation() {
    let mut truncated_viability = base_ir();
    truncated_viability.claim.objective = GameObjectiveV1::Viability;
    assert_issue(
        truncated_viability,
        GameSemanticIssueV1::StoppingSemanticsMismatch,
    );

    let mut outcome_free_stop = base_ir();
    outcome_free_stop.stopping = GameStoppingSemanticsV1::StateDependent {
        rule: GameStoppingRuleIdV1::from_bytes(bytes(54)),
    };
    assert_issue(
        outcome_free_stop,
        GameSemanticIssueV1::StoppingSemanticsMismatch,
    );

    let mut infinite_viability = base_ir();
    infinite_viability.claim.objective = GameObjectiveV1::Viability;
    infinite_viability.horizon = GameHorizonV1::Infinite {
        start: 0.0,
        unit: GameTimeUnitIdV1::from_bytes(bytes(25)),
        seconds_per_unit: 1.0,
        no_claim: GameNoClaimIdV1::from_bytes(bytes(26)),
    };
    infinite_viability.stopping = GameStoppingSemanticsV1::FirstTarget;
    assert_issue(
        infinite_viability,
        GameSemanticIssueV1::StoppingSemanticsMismatch,
    );
}

#[test]
fn hybrid_mode_information_and_zeno_scope_are_explicit() {
    let mut hybrid = base_ir();
    let modes = GameModeSetV1 {
        modes: GameModeSetIdV1::from_bytes(bytes(27)),
        model_version: hybrid.model.version,
    };
    let events = GameEventSetV1 {
        events: GameEventSetIdV1::from_bytes(bytes(28)),
        model_version: hybrid.model.version,
    };
    hybrid.model.class = GameModelClassV1::Hybrid {
        modes,
        events,
        zeno: HybridZenoScopeV1::Unresolved {
            no_claim: GameNoClaimIdV1::from_bytes(bytes(29)),
        },
    };
    hybrid.information.grants.push(InformationGrantV1 {
        player: GamePlayerV1::Controller,
        subject: InformationSubjectV1::Mode,
        observation: GameObservationMapIdV1::from_bytes(bytes(30)),
        availability: ObservationAvailabilityV1::Current,
    });
    let strategy = controller_strategy(&mut hybrid);
    strategy.representation = StrategyRepresentationV1::HybridModeFeedback {
        artifact: GameStrategyArtifactIdV1::from_bytes(bytes(31)),
    };
    strategy.dependencies.push(StrategyDependencyV1 {
        subject: InformationSubjectV1::Mode,
        access: StrategyAccessV1::Current,
    });
    let admitted = validate(hybrid.clone());
    assert_eq!(
        admitted.claim_availability(),
        GameClaimAvailabilityV1::Unknown {
            reason: GameUnknownReasonV1::ZenoUnresolved,
        }
    );

    hybrid.model.class = GameModelClassV1::Hybrid {
        modes,
        events,
        zeno: HybridZenoScopeV1::Excluded {
            witness: GameWitnessIdV1::from_bytes(bytes(32)),
        },
    };
    assert!(matches!(
        validate(hybrid).claim_availability(),
        GameClaimAvailabilityV1::Eligible { .. }
    ));

    let mut nonhybrid = base_ir();
    nonhybrid.information.grants.push(InformationGrantV1 {
        player: GamePlayerV1::Controller,
        subject: InformationSubjectV1::Mode,
        observation: GameObservationMapIdV1::from_bytes(bytes(30)),
        availability: ObservationAvailabilityV1::Current,
    });
    assert_issue(
        nonhybrid,
        GameSemanticIssueV1::HybridInformationForNonHybridModel,
    );
}

#[test]
fn set_units_model_versions_and_quantifier_domains_fail_closed() {
    let mut units = base_ir();
    units.target.context.units = GameUnitSystemIdV1::from_bytes(bytes(33));
    assert_issue(
        units,
        GameSemanticIssueV1::ContextMismatch {
            role: GameContextRoleV1::Target,
        },
    );

    let mut version = base_ir();
    version.controls.model_version = GameModelVersionIdV1::from_bytes(bytes(34));
    assert_issue(
        version,
        GameSemanticIssueV1::ContextMismatch {
            role: GameContextRoleV1::Control,
        },
    );

    let mut domain = base_ir();
    let clauses = [
        domain.quantifiers.clauses()[0],
        domain.quantifiers.clauses()[1],
        domain.quantifiers.clauses()[2],
    ];
    domain.quantifiers = GameQuantifierPrefixV1::new(vec![
        GameQuantifierClauseV1 {
            domain: GameQuantifierDomainV1::Disturbance(domain.disturbances.set),
            ..clauses[0]
        },
        clauses[1],
        clauses[2],
    ]);
    assert_issue(
        domain,
        GameSemanticIssueV1::QuantifierDomainMismatch {
            variable: GameVariableV1::Control,
        },
    );
}

#[test]
fn anticipative_policies_and_open_loop_online_reads_refuse() {
    let mut future = base_ir();
    set_disturbance_grant(&mut future, ObservationAvailabilityV1::FutureTrajectory);
    assert_issue(future, GameSemanticIssueV1::AnticipativeInformation);

    let mut policy = base_ir();
    let strategy = controller_strategy(&mut policy);
    strategy.representation = StrategyRepresentationV1::NonanticipativeFeedback {
        artifact: GameStrategyArtifactIdV1::from_bytes(bytes(35)),
        memory_states: 1,
    };
    strategy.dependencies[0].access = StrategyAccessV1::FutureTrajectory;
    assert_issue(policy, GameSemanticIssueV1::AnticipativeStrategy);

    let mut open_loop = base_ir();
    controller_strategy(&mut open_loop).representation = StrategyRepresentationV1::OpenLoop {
        artifact: GameStrategyArtifactIdV1::from_bytes(bytes(36)),
    };
    assert_issue(open_loop, GameSemanticIssueV1::OpenLoopHasOnlineDependency);
}

#[test]
fn dae_hybrid_stopping_and_unknown_strategies_are_not_overclaimed() {
    let mut dae = base_ir();
    dae.model.class = GameModelClassV1::AdmittedDae {
        index: 0,
        constraint: GameDaeConstraintIdV1::from_bytes(bytes(37)),
    };
    assert_issue(
        dae,
        GameSemanticIssueV1::InvalidValue {
            field: GameFieldV1::DaeIndex,
        },
    );

    let mut wrong_stop = base_ir();
    wrong_stop.stopping = GameStoppingSemanticsV1::HybridTerminalEvent {
        events: GameEventSetIdV1::from_bytes(bytes(38)),
    };
    assert_issue(wrong_stop, GameSemanticIssueV1::StoppingSemanticsMismatch);

    let mut unknown = base_ir();
    controller_strategy(&mut unknown).representation = StrategyRepresentationV1::Unknown {
        no_claim: GameNoClaimIdV1::from_bytes(bytes(39)),
    };
    assert_eq!(
        validate(unknown).claim_availability(),
        GameClaimAvailabilityV1::Unknown {
            reason: GameUnknownReasonV1::StrategyUnresolved,
        }
    );
}

fn parsed_problem_id(seed: u8) -> GameProblemIdV1 {
    GameProblemIdV1::parse_hex(&format!("{seed:02x}").repeat(32)).expect("typed problem id")
}

fn component(seed: u8) -> GameComponentRefV1 {
    GameComponentRefV1 {
        problem: parsed_problem_id(seed),
        model_version: GameModelVersionIdV1::from_bytes(bytes(2)),
        state_space: GameStateSpaceIdV1::from_bytes(bytes(3)),
        frame: GameFrameIdV1::from_bytes(bytes(4)),
        units: GameUnitSystemIdV1::from_bytes(bytes(5)),
    }
}

#[test]
fn parallel_composition_is_canonical_and_context_checked() {
    let mut first = base_ir();
    first.composition = GameCompositionV1::Parallel {
        components: vec![component(41), component(40)],
        interface: GameCompositionInterfaceIdV1::from_bytes(bytes(42)),
    };
    let first = validate(first);
    let mut second = base_ir();
    second.composition = GameCompositionV1::Parallel {
        components: vec![component(40), component(41)],
        interface: GameCompositionInterfaceIdV1::from_bytes(bytes(42)),
    };
    let second = validate(second);
    assert_eq!(first.identity_receipt(), second.identity_receipt());

    let mut duplicate = base_ir();
    duplicate.composition = GameCompositionV1::Parallel {
        components: vec![component(40), component(40)],
        interface: GameCompositionInterfaceIdV1::from_bytes(bytes(42)),
    };
    assert_issue(
        duplicate,
        GameSemanticIssueV1::DuplicateCompositionComponent,
    );

    let mut mismatch = base_ir();
    let mut bad = component(43);
    bad.units = GameUnitSystemIdV1::from_bytes(bytes(44));
    mismatch.composition = GameCompositionV1::Sequential {
        components: vec![component(40), bad],
        interface: GameCompositionInterfaceIdV1::from_bytes(bytes(45)),
    };
    assert_issue(
        mismatch,
        GameSemanticIssueV1::ContextMismatch {
            role: GameContextRoleV1::Component,
        },
    );
}

#[test]
fn sequential_composition_preserves_order_and_multiplicity() {
    let mut repeated_first = base_ir();
    repeated_first.composition = GameCompositionV1::Sequential {
        components: vec![component(40), component(40), component(41)],
        interface: GameCompositionInterfaceIdV1::from_bytes(bytes(55)),
    };
    let repeated_first = validate(repeated_first);

    let mut repeated_last = base_ir();
    repeated_last.composition = GameCompositionV1::Sequential {
        components: vec![component(40), component(41), component(40)],
        interface: GameCompositionInterfaceIdV1::from_bytes(bytes(55)),
    };
    let repeated_last = validate(repeated_last);

    assert_ne!(repeated_first.problem_id(), repeated_last.problem_id());
}

fn with_schema(ir: GameProblemIrV1, schema_version: u32) -> GameProblemIrV1 {
    GameProblemIrV1::with_schema_version(
        schema_version,
        ir.model,
        ir.initial,
        ir.target,
        ir.unsafe_set,
        ir.controls,
        ir.disturbances,
        ir.parameters,
        ir.claim,
        ir.quantifiers,
        ir.information,
        ir.strategies,
        ir.horizon,
        ir.stopping,
        ir.composition,
        ir.budget,
    )
}

#[test]
fn schema_caps_and_cancellation_refuse_before_publication() {
    assert_issue(
        with_schema(base_ir(), GAME_PROBLEM_SCHEMA_VERSION_V1 + 1),
        GameSemanticIssueV1::UnsupportedSchemaVersion {
            found: GAME_PROBLEM_SCHEMA_VERSION_V1 + 1,
            supported: GAME_PROBLEM_SCHEMA_VERSION_V1,
        },
    );

    let mut capped = base_ir();
    let clause = capped.quantifiers.clauses()[0];
    capped.quantifiers = GameQuantifierPrefixV1::new(vec![clause; MAX_GAME_QUANTIFIERS_V1 + 1]);
    assert_issue(
        capped,
        GameSemanticIssueV1::TooMany {
            collection: GameCollectionV1::Quantifiers,
            found: MAX_GAME_QUANTIFIERS_V1 + 1,
            limit: MAX_GAME_QUANTIFIERS_V1,
        },
    );

    let mut dependency_capped = base_ir();
    controller_strategy(&mut dependency_capped).dependencies = vec![
        StrategyDependencyV1 {
            subject: InformationSubjectV1::State,
            access: StrategyAccessV1::Current,
        };
        MAX_STRATEGY_DEPENDENCIES_V1 + 1
    ];
    assert_issue(
        dependency_capped,
        GameSemanticIssueV1::TooMany {
            collection: GameCollectionV1::StrategyDependencies,
            found: MAX_STRATEGY_DEPENDENCIES_V1 + 1,
            limit: MAX_STRATEGY_DEPENDENCIES_V1,
        },
    );

    let report = with_cx(true, |cx| {
        validate_game_problem_v1(base_ir(), cx).expect_err("pre-cancelled scope")
    });
    assert_eq!(report.issues(), &[GameSemanticIssueV1::Cancelled]);
}
