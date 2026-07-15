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
    FinalizationReceipt, HeadToHeadCharter, IdempotencyKey, LaneCharter, LaneError, MechanismId,
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

    let h2h = HeadToHeadCharter::new(lane.lane_id(), &[c1, c2], envelope(50), evidence("prereg"))
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
        HeadToHeadCharter::new(lane.lane_id(), &[c1], envelope(1), evidence("x")),
        Err(LaneError::ComparisonCandidatesInvalid)
    ));
    assert!(matches!(
        HeadToHeadCharter::new(lane.lane_id(), &[c1, c1], envelope(1), evidence("x")),
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
    let late = HeadToHeadCharter::new(lane2.lane_id(), &[d1, d2], envelope(50), evidence("late"))
        .expect("charter");
    assert!(matches!(
        fresh.preregister_comparison(late, key("late")),
        Err(LaneError::ComparisonAfterAdmission { .. })
    ));
}

/// lane-003b (comparison envelope): the declared shared budget refuses
/// the reservation that would exceed it, naming the axis.
#[test]
fn lane_003b_comparison_envelope_binds() {
    let mut ledger = PortfolioLedger::new(policy());
    let lane = charter("claim HB", "class-hb");
    let c1 = lane.mechanism_id("one", 1).expect("id");
    let c2 = lane.mechanism_id("two", 1).expect("id");
    let h2h = HeadToHeadCharter::new(lane.lane_id(), &[c1, c2], envelope(30), evidence("p"))
        .expect("charter");
    ledger
        .preregister_comparison(h2h, key("p"))
        .expect("prereg");
    ledger
        .admit(&lane, c1, envelope(20), key("c1"))
        .expect("fits");
    let refusal = ledger
        .admit(&lane, c2, envelope(20), key("c2"))
        .expect_err("20 + 20 > 30 shared budget");
    assert!(
        matches!(
            refusal,
            LaneError::ComparisonEnvelopeExceeded {
                axis: "work",
                requested: 20,
                remaining: 10
            }
        ),
        "axis-precise refusal: {refusal}"
    );
}

/// lane-004 (G0 global caps): the portfolio mechanism cap and each
/// global envelope axis bind ACROSS lanes — partitioning cannot evade
/// portfolio limits.
#[test]
fn lane_004_global_envelopes_bind_across_lanes() {
    let mut ledger = PortfolioLedger::new(PortfolioPolicy {
        global: ResourceEnvelope {
            work_units: 25,
            memory_bytes: u64::MAX,
            reviewer_slots: 100,
            falsification_capacity: 100,
        },
        max_active_mechanisms: 2,
    });
    let l1 = charter("claim 1", "class-1");
    let l2 = charter("claim 2", "class-2");
    let l3 = charter("claim 3", "class-3");
    let m1 = l1.mechanism_id("m", 1).expect("id");
    let m2 = l2.mechanism_id("m", 1).expect("id");
    let m3 = l3.mechanism_id("m", 1).expect("id");

    ledger
        .admit(
            &l1,
            m1,
            ResourceEnvelope {
                work_units: 10,
                ..Default::default()
            },
            key("m1"),
        )
        .expect("first");
    ledger
        .admit(
            &l2,
            m2,
            ResourceEnvelope {
                work_units: 10,
                ..Default::default()
            },
            key("m2"),
        )
        .expect("second");
    let capped = ledger
        .admit(
            &l3,
            m3,
            ResourceEnvelope {
                work_units: 1,
                ..Default::default()
            },
            key("m3"),
        )
        .expect_err("mechanism cap 2 binds");
    assert!(matches!(
        capped,
        LaneError::PortfolioCapExceeded { active: 2, cap: 2 }
    ));

    // Release one, then the work envelope (25) binds: 10 used + 20 > 25.
    ledger
        .finalize(
            &FinalizationReceipt::new(m2, TerminalKind::Refuted, None, evidence("r"))
                .expect("receipt"),
            key("r-m2"),
        )
        .expect("refutation releases");
    let enveloped = ledger
        .admit(
            &l3,
            m3,
            ResourceEnvelope {
                work_units: 20,
                ..Default::default()
            },
            key("m3-b"),
        )
        .expect_err("work envelope binds");
    assert!(
        matches!(
            enveloped,
            LaneError::EnvelopeExceeded {
                axis: "work",
                requested: 20,
                remaining: 15
            }
        ),
        "axis-precise: {enveloped}"
    );
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
        conflict,
        LaneError::IdempotencyConflict { original_seq: 0 }
    ));

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
    let h2h = HeadToHeadCharter::new(
        cosmetic_split.lane_id(),
        &[c1, c2],
        envelope(50),
        evidence("e"),
    )
    .expect("charter");
    ledger
        .preregister_comparison(h2h, key("h"))
        .expect("preregistration itself is fine");
    assert!(matches!(
        ledger.admit(&cosmetic_split, c1, envelope(5), key("c1")),
        Err(LaneError::IndependenceClassOccupied { .. })
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
    let bounded = first.decisions_json(2);
    assert!(
        bounded.starts_with("{\"skipped\":3,"),
        "explicit truncation: {bounded}"
    );
    // Log rows carry the fields the bead demands.
    let full = first.decisions_json(usize::MAX);
    for needle in [
        "policy_version",
        "lane",
        "mechanism",
        "idempotency",
        "request_digest",
        "verdict",
        "remedy",
    ] {
        assert!(full.contains(needle), "log field `{needle}` present");
    }
}
