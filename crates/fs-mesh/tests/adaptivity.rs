//! G0/G3 accounting tests for goal-oriented mesh evolution.

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_mesh::{
    AdaptivityAction, AdaptivityEffects, AdaptivityError, AdaptivityReceipt,
    AdaptivityReceiptAuthority, AdaptivityTrigger, BalanceStatus, ConservativeRemapError,
    LineageRecordId, MAX_REMAP_SOURCE_COVERAGE_TOLERANCE, MAX_REMAP_TARGET_CELLS, MeshStateId,
    QoiBoundSnapshot, QoiBoundTrend, QoiEvidenceId, QoiId, RemapAccounting, RemapContribution,
    RemapEvidenceId, RemapInvariantId, TopologyLineage, conservative_cell_remap,
};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x5e_ed,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn bytes(tag: u8) -> [u8; 32] {
    [tag; 32]
}

fn qoi(tag: u8) -> QoiId {
    QoiId::from_bytes(bytes(tag))
}

fn evidence(tag: u8) -> QoiEvidenceId {
    QoiEvidenceId::from_bytes(bytes(tag))
}

fn mesh(tag: u8) -> MeshStateId {
    MeshStateId::from_bytes(bytes(tag))
}

fn lineage(tag: u8) -> LineageRecordId {
    LineageRecordId::from_bytes(bytes(tag))
}

fn effects(connectivity: bool, physical_topology: bool) -> AdaptivityEffects {
    AdaptivityEffects::new(connectivity, physical_topology, true)
        .expect("fixture effects are internally consistent")
}

fn remap(
    invariant_tag: u8,
    evidence_tag: u8,
    defect: f64,
    tolerance: f64,
    projection: f64,
) -> RemapAccounting {
    RemapAccounting::new(
        RemapInvariantId::from_bytes(bytes(invariant_tag)),
        RemapEvidenceId::from_bytes(bytes(evidence_tag)),
        defect,
        tolerance,
        projection,
    )
    .expect("fixture remap accounting is admissible")
}

fn snapshot(
    state: MeshStateId,
    qoi_id: QoiId,
    evidence_tag: u8,
    estimator: f64,
    conversion: f64,
) -> QoiBoundSnapshot {
    QoiBoundSnapshot::new(state, qoi_id, evidence(evidence_tag), estimator, conversion)
        .expect("fixture accounting is admissible")
}

#[test]
fn g0_qoi_regression_is_visible_and_binds_complete_accounting() {
    let qoi_id = qoi(1);
    let transition = TopologyLineage::new(
        AdaptivityAction::AnisotropicRemesh,
        effects(true, false),
        mesh(2),
        mesh(3),
        lineage(4),
    )
    .unwrap();
    let before = snapshot(mesh(2), qoi_id, 5, 0.20, 0.05);
    let after = snapshot(mesh(3), qoi_id, 6, 0.23, 0.04);
    let remap = remap(7, 8, -2.0e-6, 1.0e-6, 3.0e-5);

    let receipt = AdaptivityReceipt::admit(
        AdaptivityTrigger::GoalOriented,
        transition,
        before,
        after,
        remap,
    )
    .unwrap();

    assert_eq!(receipt.qoi_trend(), QoiBoundTrend::Increased);
    assert!(!receipt.qoi_bound_decreased());
    assert_eq!(
        receipt.remap().balance_status(),
        BalanceStatus::ExceededDeclaredTolerance
    );
    assert_eq!(
        receipt.authority(),
        AdaptivityReceiptAuthority::DeclarationOnly
    );
    assert!(receipt.before().total_upper_bound() >= 0.25);
    assert!(receipt.after().total_upper_bound() >= 0.27);

    let json = receipt.to_json();
    assert!(json.contains("\"qoi_bound_trend\":\"increased\""));
    assert!(json.contains("\"qoi_bound_decreased\":false"));
    assert!(json.contains("\"balance_status\":\"exceeded-declared-tolerance\""));
    assert!(json.contains("\"declared_connectivity_changed\":true"));
    assert!(json.contains("\"declared_physical_topology_changed\":false"));
    assert!(json.contains("\"declared_gradient_discontinuity\":true"));
    assert!(json.contains("\"authority\":\"declaration-only\""));
    assert!(json.contains(&"04".repeat(32)));
    assert!(json.contains(&"05".repeat(32)));
    assert!(json.contains(&"06".repeat(32)));
    assert!(json.contains(&"07".repeat(32)));
    assert!(json.contains(&"08".repeat(32)));
    let expected = format!(
        concat!(
            "{{\"schema\":\"fs-mesh-adaptivity-receipt-v1\",",
            "\"authority\":\"declaration-only\",\"trigger\":\"goal-oriented\",",
            "\"action\":\"anisotropic-remesh\",\"declared_connectivity_changed\":true,",
            "\"declared_physical_topology_changed\":false,",
            "\"declared_gradient_discontinuity\":true,",
            "\"lineage_record_id\":\"{}\",\"source_mesh_state_id\":\"{}\",",
            "\"target_mesh_state_id\":\"{}\",\"qoi_id\":\"{}\",",
            "\"before_evidence_id\":\"{}\",\"before_estimator_upper_bound\":{:.17e},",
            "\"before_conversion_upper_bound\":{:.17e},\"before_total_upper_bound\":{:.17e},",
            "\"after_evidence_id\":\"{}\",\"after_estimator_upper_bound\":{:.17e},",
            "\"after_conversion_upper_bound\":{:.17e},\"after_total_upper_bound\":{:.17e},",
            "\"qoi_bound_trend\":\"increased\",\"qoi_bound_decreased\":false,",
            "\"remap_invariant_id\":\"{}\",\"remap_evidence_id\":\"{}\",",
            "\"balance_defect\":{:.17e},\"balance_tolerance\":{:.17e},",
            "\"balance_status\":\"exceeded-declared-tolerance\",",
            "\"projection_error\":{:.17e}}}"
        ),
        "04".repeat(32),
        "02".repeat(32),
        "03".repeat(32),
        "01".repeat(32),
        "05".repeat(32),
        0.20,
        0.05,
        before.total_upper_bound(),
        "06".repeat(32),
        0.23,
        0.04,
        after.total_upper_bound(),
        "07".repeat(32),
        "08".repeat(32),
        -2.0e-6,
        1.0e-6,
        3.0e-5,
    );
    assert_eq!(json, expected);
}

#[test]
fn g0_strict_decrease_and_unchanged_are_not_conflated() {
    let qoi_id = qoi(11);
    let improved = AdaptivityReceipt::admit(
        AdaptivityTrigger::Contact,
        TopologyLineage::new(
            AdaptivityAction::HRefine,
            effects(true, false),
            mesh(12),
            mesh(13),
            lineage(14),
        )
        .unwrap(),
        snapshot(mesh(12), qoi_id, 15, 0.5, 0.1),
        snapshot(mesh(13), qoi_id, 16, 0.4, 0.1),
        remap(22, 23, 0.0, 0.0, 0.0),
    )
    .unwrap();
    assert_eq!(improved.qoi_trend(), QoiBoundTrend::Decreased);
    assert!(improved.qoi_bound_decreased());
    assert_eq!(
        improved.remap().balance_status(),
        BalanceStatus::WithinDeclaredTolerance
    );

    let unchanged = AdaptivityReceipt::admit(
        AdaptivityTrigger::Wear,
        TopologyLineage::new(
            AdaptivityAction::PEnrich,
            effects(false, false),
            mesh(17),
            mesh(18),
            lineage(19),
        )
        .unwrap(),
        snapshot(mesh(17), qoi_id, 20, 0.4, 0.1),
        snapshot(mesh(18), qoi_id, 21, 0.45, 0.05),
        remap(24, 25, 1.0e-8, 1.0e-7, 2.0e-6),
    )
    .unwrap();
    assert_eq!(unchanged.qoi_trend(), QoiBoundTrend::Unchanged);
    assert!(!unchanged.qoi_bound_decreased());
}

#[test]
fn g0_qoi_trend_is_invariant_to_error_ledger_decomposition() {
    let qoi_id = qoi(26);
    let receipt = AdaptivityReceipt::admit(
        AdaptivityTrigger::GoalOriented,
        TopologyLineage::new(
            AdaptivityAction::PEnrich,
            effects(false, false),
            mesh(27),
            mesh(28),
            lineage(29),
        )
        .unwrap(),
        snapshot(mesh(27), qoi_id, 30, 1.0, 1.0),
        snapshot(mesh(28), qoi_id, 31, 2.0, 0.0),
        remap(32, 33, 0.0, 0.0, 0.0),
    )
    .unwrap();

    assert_eq!(
        receipt.before().total_upper_bound().to_bits(),
        2.0_f64.to_bits()
    );
    assert_eq!(
        receipt.after().total_upper_bound().to_bits(),
        2.0_f64.to_bits()
    );
    assert_eq!(receipt.qoi_trend(), QoiBoundTrend::Unchanged);
    assert!(!receipt.qoi_bound_decreased());
}

#[test]
fn g0_lineage_and_numerical_admission_fail_closed() {
    assert_eq!(QoiId::from_bytes([0; 32]).as_bytes(), &[0; 32]);
    assert_eq!(
        TopologyLineage::new(
            AdaptivityAction::Split,
            effects(true, true),
            mesh(31),
            mesh(31),
            lineage(32),
        ),
        Err(AdaptivityError::UnchangedMeshState)
    );
    assert!(matches!(
        QoiBoundSnapshot::new(mesh(33), qoi(34), evidence(35), f64::NAN, 0.0),
        Err(AdaptivityError::InvalidNonnegative {
            field: "estimator_upper_bound",
            ..
        })
    ));
    assert!(matches!(
        RemapAccounting::new(
            RemapInvariantId::from_bytes(bytes(42)),
            RemapEvidenceId::from_bytes(bytes(43)),
            0.0,
            0.0,
            -1.0,
        ),
        Err(AdaptivityError::InvalidNonnegative {
            field: "projection_error",
            ..
        })
    ));
    assert!(matches!(
        RemapAccounting::new(
            RemapInvariantId::from_bytes(bytes(42)),
            RemapEvidenceId::from_bytes(bytes(43)),
            0.0,
            0.0,
            f64::NAN,
        ),
        Err(AdaptivityError::InvalidNonnegative {
            field: "projection_error",
            ..
        })
    ));
    assert!(matches!(
        QoiBoundSnapshot::new(mesh(33), qoi(34), evidence(35), 0.0, -1.0),
        Err(AdaptivityError::InvalidNonnegative {
            field: "conversion_upper_bound",
            ..
        })
    ));
    assert_eq!(
        QoiBoundSnapshot::new(mesh(33), qoi(34), evidence(35), f64::MAX, f64::MAX),
        Err(AdaptivityError::QoiBoundOverflow)
    );
    assert!(matches!(
        RemapAccounting::new(
            RemapInvariantId::from_bytes(bytes(42)),
            RemapEvidenceId::from_bytes(bytes(43)),
            f64::INFINITY,
            0.0,
            0.0,
        ),
        Err(AdaptivityError::InvalidFinite {
            field: "balance_defect",
            ..
        })
    ));
    assert!(matches!(
        RemapAccounting::new(
            RemapInvariantId::from_bytes(bytes(42)),
            RemapEvidenceId::from_bytes(bytes(43)),
            0.0,
            -1.0,
            0.0,
        ),
        Err(AdaptivityError::InvalidNonnegative {
            field: "balance_tolerance",
            ..
        })
    ));

    let result = AdaptivityReceipt::admit(
        AdaptivityTrigger::Fracture,
        TopologyLineage::new(
            AdaptivityAction::Split,
            effects(true, true),
            mesh(35),
            mesh(36),
            lineage(37),
        )
        .unwrap(),
        snapshot(mesh(35), qoi(38), 39, 0.1, 0.1),
        snapshot(mesh(36), qoi(40), 41, 0.1, 0.1),
        remap(44, 45, 0.0, 0.0, 0.0),
    );
    assert_eq!(result, Err(AdaptivityError::QoiMismatch));

    let transition = TopologyLineage::new(
        AdaptivityAction::HRefine,
        effects(true, false),
        mesh(46),
        mesh(47),
        lineage(48),
    )
    .unwrap();
    let wrong_before = AdaptivityReceipt::admit(
        AdaptivityTrigger::GoalOriented,
        transition,
        snapshot(mesh(49), qoi(50), 51, 0.2, 0.1),
        snapshot(mesh(47), qoi(50), 52, 0.1, 0.1),
        remap(53, 54, 0.0, 0.0, 0.0),
    );
    assert_eq!(wrong_before, Err(AdaptivityError::BeforeStateMismatch));

    let wrong_after = AdaptivityReceipt::admit(
        AdaptivityTrigger::GoalOriented,
        transition,
        snapshot(mesh(46), qoi(50), 55, 0.2, 0.1),
        snapshot(mesh(49), qoi(50), 56, 0.1, 0.1),
        remap(57, 58, 0.0, 0.0, 0.0),
    );
    assert_eq!(wrong_after, Err(AdaptivityError::AfterStateMismatch));
}

#[test]
fn g3_receipt_bytes_replay_exactly_and_negative_zero_is_canonicalized() {
    let make = || {
        AdaptivityReceipt::admit(
            AdaptivityTrigger::MovingMesh,
            TopologyLineage::new(
                AdaptivityAction::Untangle,
                effects(false, false),
                mesh(51),
                mesh(52),
                lineage(53),
            )
            .unwrap(),
            snapshot(mesh(51), qoi(54), 55, 0.3, -0.0),
            snapshot(mesh(52), qoi(54), 56, 0.2, 0.0),
            remap(57, 58, -0.0, 0.0, -0.0),
        )
        .unwrap()
    };

    let first = make();
    let replay = make();
    assert_eq!(first, replay);
    assert_eq!(
        first.before().total_upper_bound().to_bits(),
        0.3_f64.to_bits()
    );
    assert_eq!(first.remap().balance_defect().to_bits(), 0.0_f64.to_bits());
    assert_eq!(first.to_json(), replay.to_json());
    assert!(!first.to_json().contains("-0.00000000000000000e0"));
}

#[test]
fn g0_action_discontinuity_and_domain_triggers_are_explicit() {
    assert_eq!(
        AdaptivityEffects::new(false, true, true),
        Err(AdaptivityError::PhysicalTopologyWithoutConnectivity)
    );
    assert_eq!(
        AdaptivityEffects::new(true, true, false),
        Err(AdaptivityError::GradientContinuityUnproven)
    );
    assert_eq!(
        TopologyLineage::new(
            AdaptivityAction::HRefine,
            effects(false, false),
            mesh(69),
            mesh(70),
            lineage(71),
        ),
        Err(AdaptivityError::ActionEffectsMismatch {
            action: AdaptivityAction::HRefine,
        })
    );
    assert_eq!(
        TopologyLineage::new(
            AdaptivityAction::Split,
            effects(false, false),
            mesh(69),
            mesh(70),
            lineage(71),
        ),
        Err(AdaptivityError::ActionEffectsMismatch {
            action: AdaptivityAction::Split,
        })
    );
    assert_eq!(
        TopologyLineage::new(
            AdaptivityAction::Untangle,
            effects(true, false),
            mesh(69),
            mesh(70),
            lineage(71),
        ),
        Err(AdaptivityError::ActionEffectsMismatch {
            action: AdaptivityAction::Untangle,
        })
    );

    let cases = [
        (AdaptivityAction::HRefine, "h-refine", effects(true, false)),
        (
            AdaptivityAction::HCoarsen,
            "h-coarsen",
            effects(true, false),
        ),
        (AdaptivityAction::PEnrich, "p-enrich", effects(false, false)),
        (AdaptivityAction::PReduce, "p-reduce", effects(false, false)),
        (
            AdaptivityAction::AnisotropicRemesh,
            "anisotropic-remesh",
            effects(false, false),
        ),
        (
            AdaptivityAction::Untangle,
            "untangle",
            effects(false, false),
        ),
        (AdaptivityAction::Split, "split", effects(true, false)),
        (AdaptivityAction::Merge, "merge", effects(true, true)),
    ];
    for (action, name, declared_effects) in cases {
        let transition =
            TopologyLineage::new(action, declared_effects, mesh(72), mesh(73), lineage(74))
                .unwrap();
        assert_eq!(transition.effects(), declared_effects);
        let receipt = AdaptivityReceipt::admit(
            AdaptivityTrigger::GoalOriented,
            transition,
            snapshot(mesh(72), qoi(75), 76, 0.2, 0.1),
            snapshot(mesh(73), qoi(75), 77, 0.1, 0.1),
            remap(78, 79, 0.0, 0.0, 0.0),
        )
        .unwrap();
        assert!(
            receipt
                .to_json()
                .contains(&format!("\"action\":\"{name}\""))
        );
    }

    for trigger in [
        AdaptivityTrigger::GoalOriented,
        AdaptivityTrigger::Contact,
        AdaptivityTrigger::Wear,
        AdaptivityTrigger::Fracture,
        AdaptivityTrigger::MovingMesh,
    ] {
        let receipt = AdaptivityReceipt::admit(
            trigger,
            TopologyLineage::new(
                AdaptivityAction::HRefine,
                effects(true, false),
                mesh(61),
                mesh(62),
                lineage(63),
            )
            .unwrap(),
            snapshot(mesh(61), qoi(64), 65, 0.2, 0.1),
            snapshot(mesh(62), qoi(64), 66, 0.1, 0.1),
            remap(67, 68, 0.0, 0.0, 0.0),
        )
        .unwrap();
        assert_eq!(receipt.trigger(), trigger);
    }
}

#[test]
fn g0_conservative_cell_remap_refines_extensive_values_and_feeds_accounting() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let contributions = [
            RemapContribution::new(0, 0, 0.5),
            RemapContribution::new(0, 1, 0.5),
            RemapContribution::new(1, 1, 0.25),
            RemapContribution::new(1, 2, 0.75),
        ];
        let first = conservative_cell_remap(&[2.0, 4.0], 3, &contributions, 0.0, 0.0, cx)
            .expect("binary overlap fractions conserve exactly");
        let replay = conservative_cell_remap(&[2.0, 4.0], 3, &contributions, 0.0, 0.0, cx)
            .expect("deterministic replay");
        assert_eq!(first, replay);
        assert_eq!(first.target_values(), &[1.0, 2.0, 3.0]);
        assert_eq!(first.report().source_cells(), 2);
        assert_eq!(first.report().target_cells(), 3);
        assert_eq!(first.report().contributions(), 4);
        assert_eq!(first.report().source_total(), 6.0);
        assert_eq!(first.report().target_total(), 6.0);
        assert_eq!(first.report().balance_defect().to_bits(), 0.0f64.to_bits());
        assert_eq!(first.report().maximum_source_coverage_defect(), 0.0);
        assert!(first.report().preserves_nonnegative_inputs());

        let accounting = first
            .report()
            .accounting(
                RemapInvariantId::from_bytes(bytes(80)),
                RemapEvidenceId::from_bytes(bytes(81)),
                2.5e-4,
            )
            .expect("measured result feeds the declaration-only receipt seam");
        assert_eq!(accounting.balance_defect(), 0.0);
        assert_eq!(accounting.balance_tolerance(), 0.0);
        assert_eq!(accounting.projection_error(), 2.5e-4);
        assert_eq!(
            accounting.balance_status(),
            BalanceStatus::WithinDeclaredTolerance
        );

        let json = first.report().to_json();
        assert!(json.contains("\"authority\":\"measured-f64\""));
        assert!(json.contains("\"source_cells\":2"));
        assert!(json.contains("\"target_cells\":3"));
        assert!(json.contains("\"balance_defect\":0.00000000000000000e0"));
        assert!(json.contains("caller-owned"));
    });
}

#[test]
fn g3_conservative_cell_remap_handles_signed_cancellation_without_false_positivity() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let contributions = [
            RemapContribution::new(0, 0, 0.5),
            RemapContribution::new(0, 1, 0.5),
            RemapContribution::new(1, 0, 0.5),
            RemapContribution::new(1, 1, 0.5),
        ];
        let outcome = conservative_cell_remap(&[-2.0, 2.0], 2, &contributions, 0.0, 0.0, cx)
            .expect("signed extensive quantities may cancel deterministically");
        assert_eq!(outcome.target_values(), &[0.0, 0.0]);
        assert_eq!(outcome.report().source_total(), 0.0);
        assert_eq!(outcome.report().target_total(), 0.0);
        assert!(!outcome.report().preserves_nonnegative_inputs());
    });
}

#[test]
fn g0_conservative_cell_remap_refuses_noncanonical_and_incomplete_overlap_maps() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let duplicate = [
            RemapContribution::new(0, 0, 0.5),
            RemapContribution::new(0, 0, 0.5),
        ];
        assert!(matches!(
            conservative_cell_remap(&[1.0], 1, &duplicate, 0.0, 0.0, cx),
            Err(ConservativeRemapError::NonCanonicalContributions {
                contribution: 1,
                ..
            })
        ));

        let reversed = [
            RemapContribution::new(0, 1, 0.5),
            RemapContribution::new(0, 0, 0.5),
        ];
        assert!(matches!(
            conservative_cell_remap(&[1.0], 2, &reversed, 0.0, 0.0, cx),
            Err(ConservativeRemapError::NonCanonicalContributions { .. })
        ));

        let missing_source = [
            RemapContribution::new(0, 0, 1.0),
            RemapContribution::new(2, 1, 1.0),
        ];
        assert_eq!(
            conservative_cell_remap(&[1.0, 2.0, 3.0], 2, &missing_source, 0.0, 0.0, cx),
            Err(ConservativeRemapError::MissingSourceCoverage { source: 1 })
        );

        let missing_target = [
            RemapContribution::new(0, 0, 1.0),
            RemapContribution::new(1, 0, 1.0),
        ];
        assert_eq!(
            conservative_cell_remap(&[1.0, 2.0], 2, &missing_target, 0.0, 0.0, cx),
            Err(ConservativeRemapError::MissingTargetCoverage { target: 1 })
        );

        let incomplete = [RemapContribution::new(0, 0, 0.9)];
        assert!(matches!(
            conservative_cell_remap(&[1.0], 1, &incomplete, 1.0e-12, 1.0, cx),
            Err(ConservativeRemapError::SourceCoverageExceeded { source: 0, .. })
        ));
    });
}

#[test]
fn g0_conservative_cell_remap_refuses_hostile_numbers_indices_and_tolerances() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let valid = [RemapContribution::new(0, 0, 1.0)];
        assert!(matches!(
            conservative_cell_remap(&[f64::NAN], 1, &valid, 0.0, 0.0, cx),
            Err(ConservativeRemapError::NonFiniteSourceValue { source: 0, .. })
        ));
        for fraction in [f64::NAN, -0.0, 0.0, -0.25, 1.000_000_1] {
            let invalid = [RemapContribution::new(0, 0, fraction)];
            assert!(matches!(
                conservative_cell_remap(&[1.0], 1, &invalid, 0.0, 0.0, cx),
                Err(ConservativeRemapError::InvalidFraction {
                    contribution: 0,
                    ..
                })
            ));
        }
        assert!(matches!(
            conservative_cell_remap(
                &[1.0],
                1,
                &[RemapContribution::new(1, 0, 1.0)],
                0.0,
                0.0,
                cx,
            ),
            Err(ConservativeRemapError::SourceOutOfRange { source: 1, .. })
        ));
        assert!(matches!(
            conservative_cell_remap(
                &[1.0],
                1,
                &[RemapContribution::new(0, 1, 1.0)],
                0.0,
                0.0,
                cx,
            ),
            Err(ConservativeRemapError::TargetOutOfRange { target: 1, .. })
        ));
        for tolerance in [
            f64::NAN,
            -1.0,
            MAX_REMAP_SOURCE_COVERAGE_TOLERANCE.next_up(),
        ] {
            assert!(matches!(
                conservative_cell_remap(&[1.0], 1, &valid, tolerance, 0.0, cx),
                Err(ConservativeRemapError::InvalidTolerance {
                    field: "source_coverage_tolerance",
                    ..
                })
            ));
        }
        assert!(matches!(
            conservative_cell_remap(&[1.0], 1, &valid, 0.0, f64::INFINITY, cx),
            Err(ConservativeRemapError::InvalidTolerance {
                field: "balance_tolerance",
                ..
            })
        ));
    });
}

#[test]
fn g0_conservative_cell_remap_separates_local_coverage_from_global_balance() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let fraction = 1.0 - 0.5 * MAX_REMAP_SOURCE_COVERAGE_TOLERANCE;
        let contributions = [RemapContribution::new(0, 0, fraction)];
        let error = conservative_cell_remap(
            &[1.0e9],
            1,
            &contributions,
            MAX_REMAP_SOURCE_COVERAGE_TOLERANCE,
            1.0,
            cx,
        )
        .expect_err("locally admitted coverage must not hide a large extensive imbalance");
        assert!(matches!(
            error,
            ConservativeRemapError::BalanceExceeded { .. }
        ));

        assert!(matches!(
            conservative_cell_remap(
                &[f64::MAX, f64::MAX],
                1,
                &[
                    RemapContribution::new(0, 0, 1.0),
                    RemapContribution::new(1, 0, 1.0),
                ],
                0.0,
                f64::MAX,
                cx,
            ),
            Err(ConservativeRemapError::ArithmeticOverflow {
                stage: "source-total",
                ..
            })
        ));
    });
}

#[test]
fn g4_conservative_cell_remap_cancellation_and_static_caps_refuse_before_publication() {
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    with_cx(&cancelled, |cx| {
        assert_eq!(
            conservative_cell_remap(
                &[1.0],
                1,
                &[RemapContribution::new(0, 0, 1.0)],
                0.0,
                0.0,
                cx,
            ),
            Err(ConservativeRemapError::Cancelled)
        );
    });

    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        assert_eq!(
            conservative_cell_remap(
                &[1.0],
                MAX_REMAP_TARGET_CELLS + 1,
                &[RemapContribution::new(0, 0, 1.0)],
                0.0,
                0.0,
                cx,
            ),
            Err(ConservativeRemapError::TooMany {
                field: "target_cells",
                count: MAX_REMAP_TARGET_CELLS + 1,
                cap: MAX_REMAP_TARGET_CELLS,
            })
        );
    });
}
