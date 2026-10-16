//! One-bet-per-lane admission battery (bead frankensim-ext-epic-gov-rjoq.6).
//!
//! G0 state-machine and identity laws, G3 canonicalization and
//! split/merge adversaries, exactly-once terminal release,
//! crash/retry idempotency, portfolio/comparison envelope caps, and
//! G5 deterministic decision-log replay. Each same-lane, global-cap,
//! terminal-release, and identity guard has a test that fails if the
//! guard is deleted (the in-process mutation coverage the bead
//! demands; a fault-injecting storage lane is the cross-crate E2E's
//! job and stays in the bead).

use fs_blake3::hash_domain;
use fs_govern::{
    DecisionRequest, FinalizationReceipt, HeadToHeadCharter, IdempotencyKey, LaneCharter,
    LaneError, MAX_H2H_CANDIDATES, MAX_RETAINED_DECISION_BYTES, MAX_RETAINED_DECISIONS,
    PortfolioLedger, PortfolioPolicy, ResourceEnvelope, TerminalKind,
};

fn evidence(tag: &str) -> fs_govern::ContentHash {
    hash_domain("fs-govern.test.lanes.evidence", tag.as_bytes())
}

fn charter(statement: &str, class: &str) -> LaneCharter {
    LaneCharter::new(
        statement,
        "linear elasticity, small strain, polyhedral domains",
        &["homogeneous Dirichlet boundary", "isotropic material"],
        "verified",
        "hand-checked FEEC baseline",
        "manufactured-solution refutation family",
        class,
    )
    .expect("valid charter")
}

fn envelope(work: u64) -> ResourceEnvelope {
    ResourceEnvelope {
        work_units: work,
        memory_bytes: work * 1024,
        reviewer_slots: 1,
        falsification_capacity: 1,
    }
}

fn one_axis(axis: &str, value: u64) -> ResourceEnvelope {
    let mut envelope = ResourceEnvelope::default();
    match axis {
        "work" => envelope.work_units = value,
        "memory" => envelope.memory_bytes = value,
        "reviewer" => envelope.reviewer_slots = value,
        "falsification-capacity" => envelope.falsification_capacity = value,
        _ => panic!("unknown test axis {axis}"),
    }
    envelope
}

fn one_axis_limit(axis: &str, value: u64) -> ResourceEnvelope {
    let mut limit = ResourceEnvelope {
        work_units: u64::MAX,
        memory_bytes: u64::MAX,
        reviewer_slots: u64::MAX,
        falsification_capacity: u64::MAX,
    };
    match axis {
        "work" => limit.work_units = value,
        "memory" => limit.memory_bytes = value,
        "reviewer" => limit.reviewer_slots = value,
        "falsification-capacity" => limit.falsification_capacity = value,
        _ => panic!("unknown test axis {axis}"),
    }
    limit
}

fn policy() -> PortfolioPolicy {
    PortfolioPolicy {
        global: ResourceEnvelope {
            work_units: 1_000,
            memory_bytes: 1_000 * 1024,
            reviewer_slots: 10,
            falsification_capacity: 10,
        },
        max_active_mechanisms: 8,
    }
}

/// A comparison budget that genuinely covers `candidates` reservations
/// on EVERY axis (each candidate reserves one reviewer slot and one
/// falsification run through `envelope`).
fn shared_budget(work: u64, candidates: u64) -> ResourceEnvelope {
    ResourceEnvelope {
        work_units: work,
        memory_bytes: work * 1024,
        reviewer_slots: candidates,
        falsification_capacity: candidates,
    }
}

fn key(tag: &str) -> IdempotencyKey {
    IdempotencyKey::derive(tag)
}

/// lane-001 (G3 identity): cosmetic whitespace and assumption
/// ordering/duplication collapse to ONE lane id; every semantic field
/// is mutation-sensitive.
#[test]
fn lane_001_canonical_identity() {
    let a = LaneCharter::new(
        "  div(sigma(u)) + f = 0   converges at order  p+1 ",
        "linear elasticity,   small strain, polyhedral domains",
        &[
            "isotropic material",
            "homogeneous  Dirichlet boundary",
            "isotropic material",
        ],
        "verified",
        "hand-checked FEEC baseline",
        "manufactured-solution refutation family",
        "elasticity-convergence",
    )
    .expect("charter a");
    let b = LaneCharter::new(
        "div(sigma(u)) + f = 0 converges at order p+1",
        "linear elasticity, small strain, polyhedral domains",
        &["homogeneous Dirichlet boundary", "isotropic material"],
        "verified",
        "hand-checked FEEC baseline",
        "manufactured-solution refutation family",
        "elasticity-convergence",
    )
    .expect("charter b");
    assert_eq!(a.lane_id(), b.lane_id(), "cosmetic splits collapse");
    assert_eq!(a.assumptions(), b.assumptions(), "sorted + deduped");

    // Semantic sensitivity: each field flip changes the identity.
    let base = charter("claim S", "class-x");
    let variants = [
        charter("claim S prime", "class-x"),
        LaneCharter::new(
            "claim S",
            "OTHER domain",
            &["homogeneous Dirichlet boundary", "isotropic material"],
            "verified",
            "hand-checked FEEC baseline",
            "manufactured-solution refutation family",
            "class-x",
        )
        .expect("variant"),
        LaneCharter::new(
            "claim S",
            "linear elasticity, small strain, polyhedral domains",
            &["isotropic material"],
            "verified",
            "hand-checked FEEC baseline",
            "manufactured-solution refutation family",
            "class-x",
        )
        .expect("variant"),
        LaneCharter::new(
            "claim S",
            "linear elasticity, small strain, polyhedral domains",
            &["homogeneous Dirichlet boundary", "isotropic material"],
            "estimated",
            "hand-checked FEEC baseline",
            "manufactured-solution refutation family",
            "class-x",
        )
        .expect("variant"),
        LaneCharter::new(
            "claim S",
            "linear elasticity, small strain, polyhedral domains",
            &["homogeneous Dirichlet boundary", "isotropic material"],
            "verified",
            "different baseline",
            "manufactured-solution refutation family",
            "class-x",
        )
        .expect("variant"),
        LaneCharter::new(
            "claim S",
            "linear elasticity, small strain, polyhedral domains",
            &["homogeneous Dirichlet boundary", "isotropic material"],
            "verified",
            "hand-checked FEEC baseline",
            "adversarial-mesh refutation family",
            "class-x",
        )
        .expect("variant"),
        charter("claim S", "class-y"),
    ];
    for (i, v) in variants.iter().enumerate() {
        assert_ne!(
            v.lane_id(),
            base.lane_id(),
            "field {i} must be identity-bearing"
        );
    }

    // Refusals: empty and oversized fields.
    assert!(matches!(
        LaneCharter::new("  ", "d", &[], "t", "b", "f", "c"),
        Err(LaneError::EmptyField { what: "statement" })
    ));
    let long = "x".repeat(5000);
    assert!(matches!(
        LaneCharter::new(&long, "d", &[], "t", "b", "f", "c"),
        Err(LaneError::TooLarge {
            what: "statement",
            ..
        })
    ));
}

/// lane-001b (identity authority): a mechanism remains bound to the
/// canonical lane that minted it at both comparison construction and
/// admission. A supersession cannot silently cross lanes either.
#[test]
fn lane_001b_mechanisms_are_lane_bound() {
    let lane_a = charter("claim identity A", "identity-a");
    let lane_b = charter("claim identity B", "identity-b");
    let a1 = lane_a.mechanism_id("a1", 1).expect("id");
    let a2 = lane_a.mechanism_id("a2", 1).expect("id");
    let b1 = lane_b.mechanism_id("b1", 1).expect("id");
    assert_eq!(a1.lane(), lane_a.lane_id());

    let mut ledger = PortfolioLedger::new(policy());
    let refusal = ledger
        .admit(&lane_b, a1, envelope(1), key("wrong-lane"))
        .expect_err("lane-A mechanism cannot enter lane B");
    assert!(matches!(
        refusal,
        LaneError::MechanismLaneMismatch { expected, actual }
            if expected == lane_b.lane_id() && actual == lane_a.lane_id()
    ));
    assert_eq!(ledger.active_count(), 0);

    assert!(matches!(
        HeadToHeadCharter::new(&lane_b, &[a1, a2], envelope(2), evidence("wrong-h2h")),
        Err(LaneError::MechanismLaneMismatch { .. })
    ));
    assert!(matches!(
        FinalizationReceipt::new(
            a1,
            TerminalKind::Superseded,
            Some(b1),
            evidence("cross-lane-successor")
        ),
        Err(LaneError::ReceiptInvalid { .. })
    ));

    let too_many = [a1; MAX_H2H_CANDIDATES + 1];
    assert!(matches!(
        HeadToHeadCharter::new(&lane_a, &too_many, envelope(2), evidence("too-many")),
        Err(LaneError::ComparisonCandidatesInvalid)
    ));
}

/// lane-002 (G0 same-lane guard): distinct lanes admit concurrently; a
/// second mechanism in the SAME lane refuses atomically and the ledger
/// is observably unchanged apart from the recorded refusal.
#[test]
fn lane_002_one_bet_per_lane() {
    let mut ledger = PortfolioLedger::new(policy());
    let lane_a = charter("claim A", "class-a");
    let lane_b = charter("claim B", "class-b");
    let a1 = lane_a.mechanism_id("equilibrated flux", 1).expect("id");
    let a2 = lane_a
        .mechanism_id("dual-weighted residual", 1)
        .expect("id");
    let b1 = lane_b.mechanism_id("betti witness", 1).expect("id");

    ledger
        .admit(&lane_a, a1, envelope(10), key("a1"))
        .expect("lane A admits");
    ledger
        .admit(&lane_b, b1, envelope(10), key("b1"))
        .expect("lane B admits concurrently");
    assert_eq!(ledger.active_count(), 2);

    let before_reserved = ledger.reserved();
    let refusal = ledger
        .admit(&lane_a, a2, envelope(10), key("a2"))
        .expect_err("second bet in lane A must refuse");
    assert!(
        matches!(refusal, LaneError::LaneOccupied { active, .. } if active == a1),
        "refusal names the occupant: {refusal}"
    );
    assert_eq!(ledger.active_count(), 2, "no partial admission");
    assert_eq!(ledger.reserved(), before_reserved, "no partial reservation");
    assert!(!refusal.remedy().is_empty(), "ranked remedy present");
    let last = ledger.decisions().last().expect("refusal recorded");
    assert!(!last.admitted(), "the refusal is in the log");
}

/// lane-003 (G0/G3 comparison): a preregistered head-to-head admits
/// ONLY its declared candidates under its bounded shared envelope;
/// preregistration after admission refuses; one comparison per lane.
#[test]
fn lane_003_preregistered_head_to_head() {
    let mut ledger = PortfolioLedger::new(policy());
    let lane = charter("claim H", "class-h");
    let c1 = lane.mechanism_id("candidate one", 1).expect("id");
    let c2 = lane.mechanism_id("candidate two", 1).expect("id");
    let intruder = lane.mechanism_id("undeclared intruder", 1).expect("id");

    let h2h = HeadToHeadCharter::new(&lane, &[c1, c2], shared_budget(50, 2), evidence("prereg"))
        .expect("comparison charter");
    ledger
        .preregister_comparison(h2h.clone(), key("h2h"))
        .expect("preregistration admits");
    assert!(matches!(
        ledger.preregister_comparison(h2h, key("h2h-dup")),
        Err(LaneError::ComparisonAlreadyDeclared { .. })
    ));

    ledger
        .admit(&lane, c1, envelope(20), key("c1"))
        .expect("declared candidate 1");
    ledger
        .admit(&lane, c2, envelope(20), key("c2"))
        .expect("declared candidate 2 shares the lane");
    assert!(matches!(
        ledger.admit(&lane, intruder, envelope(1), key("intruder")),
        Err(LaneError::NotADeclaredCandidate { .. })
    ));

    // Withdrawing a candidate releases its share; the terminal
    // mechanism itself can never re-enter.
    ledger
        .finalize(
            &FinalizationReceipt::new(c1, TerminalKind::Withdrawn, None, evidence("w1"))
                .expect("receipt"),
            key("w1"),
        )
        .expect("withdrawal releases");
    let refusal = ledger
        .admit(&lane, c1, envelope(20), key("c1-again"))
        .expect_err("re-admitting a terminal candidate refuses");
    assert!(matches!(refusal, LaneError::AlreadyTerminal { .. }));

    // Candidate bounds validate at construction.
    assert!(matches!(
        HeadToHeadCharter::new(&lane, &[c1], envelope(1), evidence("x")),
        Err(LaneError::ComparisonCandidatesInvalid)
    ));
    assert!(matches!(
        HeadToHeadCharter::new(&lane, &[c1, c1], envelope(1), evidence("x")),
        Err(LaneError::ComparisonCandidatesInvalid)
    ));

    // Preregistration must precede admission.
    let mut fresh = PortfolioLedger::new(policy());
    let lane2 = charter("claim H2", "class-h2");
    let d1 = lane2.mechanism_id("solo", 1).expect("id");
    let d2 = lane2.mechanism_id("late rival", 1).expect("id");
    fresh
        .admit(&lane2, d1, envelope(5), key("d1"))
        .expect("solo admits");
    let late =
        HeadToHeadCharter::new(&lane2, &[d1, d2], envelope(50), evidence("late")).expect("charter");
    assert!(matches!(
        fresh.preregister_comparison(late, key("late")),
        Err(LaneError::ComparisonAfterAdmission { .. })
    ));
}

/// lane-003b (comparison envelope): the declared shared budget refuses
/// the reservation that would exceed it, naming the axis.
#[test]
fn lane_003b_comparison_envelope_binds() {
    for axis in ["work", "memory", "reviewer", "falsification-capacity"] {
        let mut ledger = PortfolioLedger::new(policy());
        let lane = charter(
            &format!("claim comparison {axis}"),
            &format!("class-hb-{axis}"),
        );
        let c1 = lane.mechanism_id("one", 1).expect("id");
        let c2 = lane.mechanism_id("two", 1).expect("id");
        let h2h = HeadToHeadCharter::new(
            &lane,
            &[c1, c2],
            one_axis(axis, 3),
            evidence(&format!("p-{axis}")),
        )
        .expect("charter");
        ledger
            .preregister_comparison(h2h, key(&format!("p-{axis}")))
            .expect("prereg");
        ledger
            .admit(&lane, c1, one_axis(axis, 2), key(&format!("c1-{axis}")))
            .expect("first candidate fits");
        let refusal = ledger
            .admit(&lane, c2, one_axis(axis, 2), key(&format!("c2-{axis}")))
            .expect_err("2 + 2 exceeds the shared cap 3");
        assert!(
            matches!(
                refusal,
                LaneError::ComparisonEnvelopeExceeded {
                    axis: refused_axis,
                    requested: 2,
                    remaining: 1
                } if refused_axis == axis
            ),
            "axis-precise comparison refusal for {axis}: {refusal}"
        );
    }
}

/// lane-004 (G0 global caps): the portfolio mechanism cap and each
/// global envelope axis bind ACROSS lanes — partitioning cannot evade
/// portfolio limits.
#[test]
fn lane_004_global_envelopes_bind_across_lanes() {
    for axis in ["work", "memory", "reviewer", "falsification-capacity"] {
        let mut ledger = PortfolioLedger::new(PortfolioPolicy {
            global: one_axis_limit(axis, 10),
            max_active_mechanisms: 8,
        });
        let first = charter(&format!("claim first {axis}"), &format!("first-{axis}"));
        let second = charter(&format!("claim second {axis}"), &format!("second-{axis}"));
        let m1 = first.mechanism_id("m", 1).expect("id");
        let m2 = second.mechanism_id("m", 1).expect("id");
        ledger
            .admit(&first, m1, one_axis(axis, 6), key(&format!("first-{axis}")))
            .expect("first reservation fits");
        let refusal = ledger
            .admit(
                &second,
                m2,
                one_axis(axis, 5),
                key(&format!("second-{axis}")),
            )
            .expect_err("6 + 5 exceeds global cap 10");
        assert!(
            matches!(
                refusal,
                LaneError::EnvelopeExceeded {
                    axis: refused_axis,
                    requested: 5,
                    remaining: 4
                } if refused_axis == axis
            ),
            "axis-precise global refusal for {axis}: {refusal}"
        );
    }

    let mut capped = PortfolioLedger::new(PortfolioPolicy {
        global: one_axis_limit("work", u64::MAX),
        max_active_mechanisms: 1,
    });
    let first = charter("cap claim first", "cap-first");
    let second = charter("cap claim second", "cap-second");
    let m1 = first.mechanism_id("m", 1).expect("id");
    let m2 = second.mechanism_id("m", 1).expect("id");
    capped
        .admit(&first, m1, ResourceEnvelope::default(), key("cap-first"))
        .expect("first mechanism");
    assert!(matches!(
        capped.admit(&second, m2, ResourceEnvelope::default(), key("cap-second")),
        Err(LaneError::PortfolioCapExceeded { active: 1, cap: 1 })
    ));
}

/// lane-005 (G0 terminal release): a slot releases EXACTLY ONCE
/// against a valid receipt; stalled work never auto-releases; terminal
/// is permanent; zero-evidence and mismatched receipts refuse; a
/// supersession names a distinct successor.
#[test]
fn lane_005_terminal_release_exactly_once() {
    let mut ledger = PortfolioLedger::new(policy());
    let lane = charter("claim T", "class-t");
    let m = lane.mechanism_id("the bet", 1).expect("id");
    let successor = lane.mechanism_id("the bet", 2).expect("id");
    ledger
        .admit(&lane, m, envelope(10), key("m"))
        .expect("admit");

    // Stalled/Unknown never silently releases: the lane stays occupied
    // no matter how many admission attempts arrive.
    for i in 0..3 {
        let rival = lane.mechanism_id("rival", i).expect("id");
        assert!(matches!(
            ledger.admit(&lane, rival, envelope(1), key(&format!("rival{i}"))),
            Err(LaneError::LaneOccupied { .. })
        ));
    }

    // Receipt validation: zero evidence, self-supersession, missing
    // successor, spurious successor.
    assert!(matches!(
        FinalizationReceipt::new(
            m,
            TerminalKind::Refuted,
            None,
            fs_govern::ContentHash([0; 32])
        ),
        Err(LaneError::ReceiptInvalid { .. })
    ));
    assert!(matches!(
        FinalizationReceipt::new(m, TerminalKind::Superseded, Some(m), evidence("s")),
        Err(LaneError::ReceiptInvalid { .. })
    ));
    assert!(matches!(
        FinalizationReceipt::new(m, TerminalKind::Superseded, None, evidence("s")),
        Err(LaneError::ReceiptInvalid { .. })
    ));
    assert!(matches!(
        FinalizationReceipt::new(m, TerminalKind::Withdrawn, Some(successor), evidence("s")),
        Err(LaneError::ReceiptInvalid { .. })
    ));

    // Finalizing a mechanism that was never admitted refuses.
    let ghost = lane.mechanism_id("ghost", 9).expect("id");
    assert!(matches!(
        ledger.finalize(
            &FinalizationReceipt::new(ghost, TerminalKind::Withdrawn, None, evidence("g"))
                .expect("receipt"),
            key("ghost"),
        ),
        Err(LaneError::UnknownMechanism { .. })
    ));
    assert_eq!(
        ledger
            .decisions()
            .last()
            .expect("ghost refusal logged")
            .lane,
        ghost.lane(),
        "unknown finalization retains the mechanism's validated lane, not a forged placeholder"
    );

    // Valid supersession releases exactly once.
    let receipt = FinalizationReceipt::new(
        m,
        TerminalKind::Superseded,
        Some(successor),
        evidence("sup"),
    )
    .expect("receipt");
    ledger.finalize(&receipt, key("sup")).expect("releases");
    assert_eq!(ledger.active_count(), 0);
    assert_eq!(
        ledger.reserved(),
        ResourceEnvelope::default(),
        "capacity returned once"
    );

    // Replay of the SAME finalize is idempotent-Ok and does not
    // double-release; a NEW finalize on the terminal mechanism refuses.
    ledger
        .finalize(&receipt, key("sup"))
        .expect("idempotent replay");
    assert_eq!(ledger.reserved(), ResourceEnvelope::default());
    let again =
        FinalizationReceipt::new(m, TerminalKind::Withdrawn, None, evidence("w")).expect("receipt");
    assert!(matches!(
        ledger.finalize(&again, key("again")),
        Err(LaneError::AlreadyTerminal { .. })
    ));

    // Terminal is permanent: re-admission refuses; the successor may
    // now take the lane.
    assert!(matches!(
        ledger.admit(&lane, m, envelope(1), key("m-again")),
        Err(LaneError::AlreadyTerminal { .. })
    ));
    ledger
        .admit(&lane, successor, envelope(10), key("succ"))
        .expect("successor admits");
}

/// lane-006 (G4 crash/retry): idempotent replays return the recorded
/// decision without double-charging; a different request under a used
/// key refuses with the original sequence named.
#[test]
fn lane_006_idempotency() {
    let mut ledger = PortfolioLedger::new(policy());
    let lane = charter("claim I", "class-i");
    let m = lane.mechanism_id("m", 1).expect("id");
    ledger
        .admit(&lane, m, envelope(10), key("k"))
        .expect("admit");
    let after_first = (
        ledger.active_count(),
        ledger.reserved(),
        ledger.decisions().len(),
    );

    // Byte-identical retry (crash between commit and ack): same Ok, no
    // second charge, ONE decision row.
    ledger
        .admit(&lane, m, envelope(10), key("k"))
        .expect("replay is Ok");
    assert_eq!(
        (
            ledger.active_count(),
            ledger.reserved(),
            ledger.decisions().len()
        ),
        after_first,
        "replay neither charges nor re-records"
    );

    // Same key, different request: refuse naming the original.
    let other = lane.mechanism_id("other", 1).expect("id");
    let conflict = ledger
        .admit(&lane, other, envelope(10), key("k"))
        .expect_err("key reuse for a new request refuses");
    assert!(matches!(
        &conflict,
        LaneError::IdempotencyConflict { original_seq: 0 }
    ));
    assert_eq!(
        ledger.decisions().len(),
        after_first.2 + 1,
        "conflicting reuse is itself auditable"
    );
    assert!(matches!(
        ledger
            .decisions()
            .last()
            .and_then(|decision| decision.refusal.as_ref()),
        Some(LaneError::IdempotencyConflict { original_seq: 0 })
    ));
    let rows_after_conflict = ledger.decisions().len();
    let replayed_conflict = ledger
        .admit(&lane, other, envelope(10), key("k"))
        .expect_err("exact retry of the conflicting request replays");
    assert_eq!(replayed_conflict, conflict);
    assert_eq!(
        ledger.decisions().len(),
        rows_after_conflict,
        "exact conflict replay does not consume another bounded row"
    );

    // Refusals replay too: the recorded refusal returns without a new row.
    let rival = lane.mechanism_id("rival", 1).expect("id");
    let e1 = ledger
        .admit(&lane, rival, envelope(1), key("rv"))
        .expect_err("occupied");
    let rows = ledger.decisions().len();
    let e2 = ledger
        .admit(&lane, rival, envelope(1), key("rv"))
        .expect_err("replayed refusal");
    assert_eq!(e1, e2, "replay returns the recorded refusal");
    assert_eq!(ledger.decisions().len(), rows, "no duplicate decision row");
}

/// lane-007 (G3 split adversary): two textually different lanes that
/// DECLARE the same independence class share one bet; a preregistered
/// comparison on another lane cannot evade the backstop; distinct
/// classes stay independent.
#[test]
fn lane_007_independence_class_backstop() {
    let mut ledger = PortfolioLedger::new(policy());
    let original = charter("claim S", "shared-fate");
    let cosmetic_split = charter("claim S, but rephrased as a split", "shared-fate");
    let independent = charter("claim genuinely elsewhere", "other-fate");
    let m1 = original.mechanism_id("m", 1).expect("id");
    let m2 = cosmetic_split.mechanism_id("m", 1).expect("id");
    let m3 = independent.mechanism_id("m", 1).expect("id");
    assert_ne!(
        original.lane_id(),
        cosmetic_split.lane_id(),
        "different lanes..."
    );

    ledger
        .admit(&original, m1, envelope(10), key("m1"))
        .expect("first bet");
    let blocked = ledger
        .admit(&cosmetic_split, m2, envelope(10), key("m2"))
        .expect_err("...but one declared falsification fate = one bet");
    assert!(matches!(blocked, LaneError::IndependenceClassOccupied { active } if active == m1));
    ledger
        .admit(&independent, m3, envelope(10), key("m3"))
        .expect("distinct class admits");

    // Comparison evasion: preregistering a comparison on the split
    // lane does not bypass the class backstop.
    let c1 = cosmetic_split.mechanism_id("c1", 1).expect("id");
    let c2 = cosmetic_split.mechanism_id("c2", 1).expect("id");
    let h2h = HeadToHeadCharter::new(&cosmetic_split, &[c1, c2], envelope(50), evidence("e"))
        .expect("charter");
    ledger
        .preregister_comparison(h2h, key("h"))
        .expect("preregistration itself is fine");
    assert!(matches!(
        ledger.admit(&cosmetic_split, c1, envelope(5), key("c1")),
        Err(LaneError::IndependenceClassOccupied { .. })
    ));

    // Regression: every surviving candidate, not just one map
    // representative, keeps the class occupied.
    let mut survivors = PortfolioLedger::new(policy());
    let compared = charter("claim compared", "survivor-fate");
    let split = charter("claim compared but split", "survivor-fate");
    let first = compared.mechanism_id("first", 1).expect("id");
    let second = compared.mechanism_id("second", 1).expect("id");
    let split_mechanism = split.mechanism_id("split", 1).expect("id");
    let comparison = HeadToHeadCharter::new(
        &compared,
        &[first, second],
        shared_budget(30, 2),
        evidence("survivor-comparison"),
    )
    .expect("comparison");
    survivors
        .preregister_comparison(comparison, key("survivor-prereg"))
        .expect("preregister");
    survivors
        .admit(&compared, first, envelope(10), key("survivor-first"))
        .expect("first candidate");
    survivors
        .admit(&compared, second, envelope(10), key("survivor-second"))
        .expect("second candidate");
    survivors
        .finalize(
            &FinalizationReceipt::new(
                first,
                TerminalKind::Withdrawn,
                None,
                evidence("survivor-first-final"),
            )
            .expect("receipt"),
            key("survivor-first-final"),
        )
        .expect("first candidate finalizes");
    assert!(matches!(
        survivors.admit(
            &split,
            split_mechanism,
            envelope(1),
            key("survivor-split")
        ),
        Err(LaneError::IndependenceClassOccupied { active }) if active == second
    ));
}

/// lane-008 (G5 replay): the same request sequence in a fresh ledger
/// reproduces the decision log EXACTLY (JSON-identical), and the
/// bounded emitter names what it skipped.
#[test]
fn lane_008_deterministic_replay() {
    let run = || {
        let mut ledger = PortfolioLedger::new(policy());
        let lane_a = charter("claim A", "class-a");
        let lane_b = charter("claim B", "class-b");
        let a1 = lane_a.mechanism_id("m", 1).expect("id");
        let a2 = lane_a.mechanism_id("n", 1).expect("id");
        let b1 = lane_b.mechanism_id("m", 1).expect("id");
        let _ = ledger.admit(&lane_a, a1, envelope(10), key("1"));
        let _ = ledger.admit(&lane_a, a2, envelope(10), key("2")); // refused
        let _ = ledger.admit(&lane_b, b1, envelope(10), key("3"));
        let _ = ledger.finalize(
            &FinalizationReceipt::new(a1, TerminalKind::Refuted, None, evidence("r"))
                .expect("receipt"),
            key("4"),
        );
        let _ = ledger.admit(&lane_a, a2, envelope(10), key("5")); // now admits
        ledger
    };
    let first = run();
    let second = run();
    assert_eq!(first, second, "whole-ledger determinism");
    assert_eq!(
        first.decisions_json(usize::MAX),
        second.decisions_json(usize::MAX),
        "JSON log replay-identical"
    );
    assert_eq!(first.decisions().len(), 5);
    assert_eq!(first.decisions()[0].policy, policy());
    match &first.decisions()[0].request {
        DecisionRequest::Admit {
            charter,
            reservation,
        } => {
            assert_eq!(charter.statement(), "claim A");
            assert_eq!(
                charter.admissible_domain(),
                "linear elasticity, small strain, polyhedral domains"
            );
            assert_eq!(charter.assumptions().len(), 2);
            assert_eq!(*reservation, envelope(10));
        }
        other => panic!("expected replayable admit request, got {other:?}"),
    }
    match &first.decisions()[3].request {
        DecisionRequest::Finalize { kind, released, .. } => {
            assert_eq!(*kind, TerminalKind::Refuted);
            assert_eq!(*released, Some(envelope(10)));
        }
        other => panic!("expected replayable finalize request, got {other:?}"),
    }
    let bounded = first.decisions_json(2);
    assert!(
        bounded.starts_with("{\"skipped\":3,"),
        "explicit truncation: {bounded}"
    );
    // Log rows carry the fields the bead demands.
    let full = first.decisions_json(usize::MAX);
    for needle in [
        "policy_version",
        "max_active_mechanisms",
        "max_retained_decisions",
        "lane",
        "mechanism",
        "mechanism_lane",
        "idempotency",
        "request_digest",
        "statement",
        "admissible_domain",
        "assumptions",
        "reservation",
        "falsification_capacity",
        "terminal_kind",
        "released",
        "verdict",
        "remedy",
    ] {
        assert!(full.contains(needle), "log field `{needle}` present");
    }
}

/// lane-009 (bounded retention): fresh traffic refuses before either
/// retained decisions or idempotency bindings can grow without bound.
/// One slot/key remains reserved for every active mechanism's terminal
/// transition, and exact retries still replay when the ledger is full.
#[test]
fn lane_009_retention_is_bounded_and_finalization_stays_available() {
    let mut ledger = PortfolioLedger::new(policy());
    let active_lane = charter("retention active claim", "retention-active");
    let wrong_lane = charter("retention wrong claim", "retention-wrong");
    let mechanism = active_lane.mechanism_id("active", 1).expect("id");
    ledger
        .admit(
            &active_lane,
            mechanism,
            envelope(1),
            key("retention-active"),
        )
        .expect("active mechanism");

    let mut hit_cap = false;
    for index in 0..(MAX_RETAINED_DECISIONS * 2) {
        match ledger.admit(
            &wrong_lane,
            mechanism,
            ResourceEnvelope::default(),
            key(&format!("retention-wrong-{index}")),
        ) {
            Err(LaneError::MechanismLaneMismatch { .. }) => {}
            Err(LaneError::RetentionCapacityExceeded {
                axis: "decision-count",
                ..
            }) => {
                hit_cap = true;
                break;
            }
            other => panic!("unexpected retention result: {other:?}"),
        }
    }
    assert!(hit_cap, "fresh decisions eventually refuse at the hard cap");
    assert_eq!(
        ledger.decisions().len() + usize::try_from(ledger.active_count()).expect("small"),
        MAX_RETAINED_DECISIONS,
        "one decision slot is reserved for the active mechanism"
    );
    assert!(ledger.retained_decision_bytes() <= MAX_RETAINED_DECISION_BYTES);

    let receipt = FinalizationReceipt::new(
        mechanism,
        TerminalKind::Withdrawn,
        None,
        evidence("retention-final"),
    )
    .expect("receipt");
    ledger
        .finalize(&receipt, key("retention-final"))
        .expect("reserved finalization slot remains usable");
    assert_eq!(ledger.decisions().len(), MAX_RETAINED_DECISIONS);
    ledger
        .finalize(&receipt, key("retention-final"))
        .expect("exact retry replays even at capacity");
    assert!(matches!(
        ledger.admit(
            &wrong_lane,
            wrong_lane.mechanism_id("new", 1).expect("id"),
            ResourceEnvelope::default(),
            key("after-full")
        ),
        Err(LaneError::RetentionCapacityExceeded { .. })
    ));
}
