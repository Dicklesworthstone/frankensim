//! G0 battery for the cross-repository compatibility suite and release train.
//!
//! The registry tests pin structural integrity: every registered test belongs
//! to a claim its own surface declares, and every uncovered surface says why.
//! The delta tests use the REAL 2026-07-24 asupersync drift as their fixture,
//! because that incident is what this machinery exists to name. The bump tests
//! are all negative: an unrun suite, a failing suite, an empty execution, a
//! missing result, and an uncovered moved sibling must each refuse, and a
//! refusal reports every reason at once rather than one per attempt.

use fs_govern::compatibility::{
    BumpAttempt, BumpKind, BumpRefusal, BumpVerdict, EmergencyClass, GoldenDisposition, PinDelta,
    PinMovement, PinRow, ReviewPriority, SURFACES, SuiteOutcome, SuiteResult,
    coupled_golden_surfaces, evaluate_bump, moved, pin_delta, render_registry, surface,
};

/// The recorded lock as of 2026-07-24.
fn recorded() -> Vec<PinRow> {
    vec![
        PinRow::new(
            "asupersync",
            "0.3.8",
            "5973b0ff31f405ae90fa9e6e2ef5f61a75c5b78b",
        ),
        PinRow::new(
            "frankensqlite",
            "unversioned-workspace",
            "987cfb97f86d3fca4d9b44e7871f427636b10126",
        ),
    ]
}

/// Green results for every surface that is required on EVERY train: the
/// P1/P2 siblings that sit in the runtime graph.
fn fully_evidenced() -> Vec<SuiteResult> {
    vec![
        green("asupersync"),
        green("franken_numpy"),
        green("frankensqlite"),
    ]
}

fn green(lib: &str) -> SuiteResult {
    SuiteResult::new(
        lib,
        SuiteOutcome::Executed {
            passed: 6,
            failed: 0,
        },
    )
}

/// Every registered surface is internally consistent: a test may only cite a
/// claim its own surface declares, and a surface without tests must say why.
#[test]
fn registry_is_structurally_consistent() {
    assert!(!SURFACES.is_empty());

    let mut libs: Vec<&str> = SURFACES.iter().map(|entry| entry.lib).collect();
    let sorted = {
        let mut copy = libs.clone();
        copy.sort_unstable();
        copy
    };
    assert_eq!(libs, sorted, "surfaces must be in canonical lib order");
    libs.dedup();
    assert_eq!(libs.len(), SURFACES.len(), "no duplicate sibling");

    for entry in SURFACES {
        for test in entry.tests {
            assert!(
                entry.claims.contains(&test.claim),
                "{} test {} cites claim {:?} that its surface does not declare",
                entry.lib,
                test.test_name,
                test.claim
            );
            assert!(!test.crate_name.is_empty());
            assert!(!test.test_target.is_empty());
            assert!(!test.test_name.is_empty());
        }
        if entry.tests.is_empty() {
            assert!(
                entry.no_test_reason.is_some(),
                "{} has no tests and no recorded reason; an unexplained gap is not admissible",
                entry.lib
            );
            assert!(
                entry.selectors().is_empty(),
                "{} must not render a runnable selector with no tests",
                entry.lib
            );
        } else {
            assert!(
                entry.no_test_reason.is_none(),
                "{} has tests, so a no-test reason is contradictory",
                entry.lib
            );
            let selectors = entry.selectors();
            assert!(!selectors.is_empty(), "{} renders a selector", entry.lib);
            for selector in &selectors {
                assert!(selector.starts_with("cargo test --locked -p "));
                // A unit-test surface selects with --lib; an integration one
                // with --test. Every group must do one or the other, or the
                // invocation silently runs the whole crate.
                assert!(
                    selector.contains(" --lib ") || selector.contains(" --test "),
                    "{} selector runs the whole crate: {selector}",
                    entry.lib
                );
                // The invocation must run EXACTLY the registered claim set, or
                // a test outside it can fail the surface for a reason the
                // registry never claimed.
                assert!(
                    selector.contains(" -- --exact "),
                    "{} selector does not pin its test set: {selector}",
                    entry.lib
                );
            }
            // A feature-gated boundary MUST carry --features, or the selector
            // compiles the surface out and reports a vacuous pass.
            for test in entry.tests {
                if !test.required_features.is_empty() {
                    assert!(
                        selectors.iter().any(|selector| test
                            .required_features
                            .iter()
                            .all(|feature| selector.contains(feature))),
                        "{} test {} needs features {:?} that no selector enables",
                        entry.lib,
                        test.test_name,
                        test.required_features
                    );
                }
            }
        }
    }

    // The two critical siblings are the ones that must carry real coverage.
    for lib in ["asupersync", "frankensqlite"] {
        let entry = surface(lib).expect("critical sibling is registered");
        assert_eq!(entry.priority, ReviewPriority::P1);
        assert!(
            entry.test_count() >= 4,
            "{lib} carries only {} tests",
            entry.test_count()
        );
    }

    // Exactly one surface is uncovered, and it is the pinned-unused one.
    let uncovered: Vec<&str> = SURFACES
        .iter()
        .filter(|entry| entry.tests.is_empty())
        .map(|entry| entry.lib)
        .collect();
    assert_eq!(
        uncovered,
        vec!["frankenpandas"],
        "every sibling with a consumer must carry boundary coverage"
    );

    // A pinned-unused sibling makes no claim at all.
    let pandas = surface("frankenpandas").expect("registered");
    assert_eq!(pandas.priority, ReviewPriority::P4);
    assert!(pandas.claims.is_empty());
    assert!(!pandas.priority.required_every_train());
}

/// The 2026-07-24 incident: the aggregate hash said only "a constellation repo
/// moved". The delta must NAME asupersync and show both endpoints.
#[test]
fn pin_delta_names_the_sibling_that_moved() {
    let live = vec![
        PinRow::new(
            "asupersync",
            "0.3.9",
            "2691a88080d74b353405427d9af48d10e077699d",
        ),
        PinRow::new(
            "frankensqlite",
            "unversioned-workspace",
            "987cfb97f86d3fca4d9b44e7871f427636b10126",
        ),
    ];
    let deltas = pin_delta(&recorded(), &live);
    assert_eq!(
        deltas.len(),
        2,
        "every sibling is reported, not just movers"
    );

    let movers = moved(&deltas);
    assert_eq!(movers.len(), 1);
    assert_eq!(movers[0].lib, "asupersync");

    let description = movers[0].describe();
    assert!(description.contains("asupersync MOVED"), "{description}");
    assert!(description.contains("0.3.8@5973b0ff"), "{description}");
    assert!(description.contains("0.3.9@2691a880"), "{description}");

    // The sibling that did NOT move is reported as unchanged, not omitted.
    let quiet = deltas
        .iter()
        .find(|delta| delta.lib == "frankensqlite")
        .expect("present");
    assert_eq!(quiet.movement, PinMovement::Unchanged);
    assert!(!quiet.movement.is_movement());
}

/// A sibling appearing or disappearing is reported, never silently ignored.
#[test]
fn pin_delta_reports_added_and_removed_siblings() {
    let candidate = vec![
        PinRow::new(
            "asupersync",
            "0.3.8",
            "5973b0ff31f405ae90fa9e6e2ef5f61a75c5b78b",
        ),
        PinRow::new("newcomer", "0.1.0", "abc"),
    ];
    let deltas = pin_delta(&recorded(), &candidate);
    let by_lib = |lib: &str| {
        deltas
            .iter()
            .find(|delta| delta.lib == lib)
            .map(|delta| delta.movement.clone())
            .expect("present")
    };
    assert_eq!(by_lib("newcomer"), PinMovement::Added);
    assert_eq!(by_lib("frankensqlite"), PinMovement::Removed);
    assert_eq!(by_lib("asupersync"), PinMovement::Unchanged);
    assert_eq!(moved(&deltas).len(), 2);
}

/// An unrun suite is never a pass, and neither is an execution of zero tests.
#[test]
fn unrun_and_empty_suites_are_never_green() {
    assert!(!SuiteOutcome::NotRun.is_green());
    assert!(
        !SuiteOutcome::Executed {
            passed: 0,
            failed: 0
        }
        .is_green(),
        "a selector that matched nothing is not evidence"
    );
    assert!(
        !SuiteOutcome::Executed {
            passed: 5,
            failed: 1
        }
        .is_green()
    );
    assert!(
        SuiteOutcome::Executed {
            passed: 1,
            failed: 0
        }
        .is_green()
    );
}

fn attempt_with(results: Vec<SuiteResult>) -> BumpAttempt {
    BumpAttempt {
        kind: BumpKind::ReleaseTrain,
        deltas: vec![PinDelta {
            lib: "asupersync".to_string(),
            movement: PinMovement::Moved {
                from_version: "0.3.8".to_string(),
                from_head: "5973b0ff".to_string(),
                to_version: "0.3.9".to_string(),
                to_head: "2691a880".to_string(),
            },
        }],
        results,
        golden: GoldenDisposition::NoGoldenSurface,
        coupled_goldens: Vec::new(),
    }
}

/// The core gate: a bump cannot land without executed, green evidence.
#[test]
fn bump_is_refused_without_green_evidence() {
    // Not run at all.
    let verdict = evaluate_bump(&attempt_with(vec![
        SuiteResult::new("asupersync", SuiteOutcome::NotRun),
        green("frankensqlite"),
    ]));
    assert!(!verdict.admitted());
    let BumpVerdict::Refused { reasons } = verdict else {
        panic!("expected refusal")
    };
    assert!(reasons.contains(&BumpRefusal::NotExecuted {
        lib: "asupersync".to_string()
    }));

    // Executed with failures.
    let verdict = evaluate_bump(&attempt_with(vec![
        SuiteResult::new(
            "asupersync",
            SuiteOutcome::Executed {
                passed: 5,
                failed: 2,
            },
        ),
        green("frankensqlite"),
    ]));
    let BumpVerdict::Refused { reasons } = verdict else {
        panic!("expected refusal")
    };
    assert!(reasons.contains(&BumpRefusal::FailingSurface {
        lib: "asupersync".to_string(),
        failed: 2
    }));

    // Executed but empty.
    let verdict = evaluate_bump(&attempt_with(vec![
        SuiteResult::new(
            "asupersync",
            SuiteOutcome::Executed {
                passed: 0,
                failed: 0,
            },
        ),
        green("frankensqlite"),
    ]));
    let BumpVerdict::Refused { reasons } = verdict else {
        panic!("expected refusal")
    };
    assert!(reasons.contains(&BumpRefusal::EmptyExecution {
        lib: "asupersync".to_string()
    }));

    // No result reported at all.
    let verdict = evaluate_bump(&attempt_with(vec![green("frankensqlite")]));
    let BumpVerdict::Refused { reasons } = verdict else {
        panic!("expected refusal")
    };
    assert!(reasons.contains(&BumpRefusal::MissingResult {
        lib: "asupersync".to_string()
    }));
}

/// A sibling that did not move still blocks the train when it is P1/P2:
/// one sibling's move can break another's surface.
#[test]
fn unmoved_critical_surfaces_are_still_required() {
    // frankensqlite did not move, but its evidence is missing.
    let verdict = evaluate_bump(&attempt_with(vec![green("asupersync")]));
    let BumpVerdict::Refused { reasons } = verdict else {
        panic!("expected refusal")
    };
    assert!(reasons.contains(&BumpRefusal::MissingResult {
        lib: "frankensqlite".to_string()
    }));

    // A P3 sibling that did not move is not required.
    assert!(!reasons.iter().any(|reason| matches!(
        reason,
        BumpRefusal::MissingResult { lib } if lib == "frankentorch"
    )));
}

/// A moved sibling that ASSERTS NOTHING needs no evidence; a moved sibling
/// nobody registered cannot be adjudicated at all.
#[test]
fn vacuous_movers_are_admitted_and_unregistered_ones_are_refused() {
    // frankenpandas is pinned-unused: no claims, no runtime consumer. A pin
    // move cannot break something nothing depends on, so demanding evidence
    // would be demanding proof of a vacuous claim.
    let pandas = surface("frankenpandas").expect("registered");
    assert!(pandas.claims.is_empty() && pandas.runtime_consumers.is_empty());

    let mut vacuous = attempt_with(fully_evidenced());
    vacuous.deltas.push(PinDelta {
        lib: "frankenpandas".to_string(),
        movement: PinMovement::Moved {
            from_version: "0.1.2".to_string(),
            from_head: "803efc1c".to_string(),
            to_version: "0.1.2".to_string(),
            to_head: "2dded976".to_string(),
        },
    });
    assert!(
        evaluate_bump(&vacuous).admitted(),
        "a sibling asserting nothing cannot require evidence"
    );

    // But a sibling nobody registered is refused: it cannot be adjudicated,
    // and silence is not the same as "asserts nothing".
    let mut unknown = attempt_with(fully_evidenced());
    unknown.deltas.push(PinDelta {
        lib: "mystery-lib".to_string(),
        movement: PinMovement::Added,
    });
    let BumpVerdict::Refused { reasons } = evaluate_bump(&unknown) else {
        panic!("expected refusal")
    };
    assert!(reasons.contains(&BumpRefusal::UnregisteredSibling {
        lib: "mystery-lib".to_string()
    }));
}

/// The exemption above must not become a laundering hole: it requires BOTH no
/// claims AND no runtime consumer, and today every sibling that declares a
/// claim carries tests for it.
#[test]
fn every_claiming_sibling_carries_coverage() {
    for entry in SURFACES {
        if entry.claims.is_empty() {
            assert!(
                entry.runtime_consumers.is_empty(),
                "{} declares no claim yet is consumed at runtime — it cannot be exempt",
                entry.lib
            );
            assert!(entry.tests.is_empty());
        } else {
            assert!(
                !entry.tests.is_empty(),
                "{} declares claims but registers no test for them",
                entry.lib
            );
        }
    }
}

/// A bump with nothing to bump is refused, so an empty train cannot be
/// recorded as a successful one.
#[test]
fn a_bump_with_no_movement_is_refused() {
    let attempt = BumpAttempt {
        kind: BumpKind::ReleaseTrain,
        deltas: vec![PinDelta {
            lib: "asupersync".to_string(),
            movement: PinMovement::Unchanged,
        }],
        results: fully_evidenced(),
        golden: GoldenDisposition::NoGoldenSurface,
        coupled_goldens: Vec::new(),
    };
    let BumpVerdict::Refused { reasons } = evaluate_bump(&attempt) else {
        panic!("expected refusal")
    };
    assert!(reasons.contains(&BumpRefusal::NoMovement));
}

/// The admitted path, and the fact that refusal reports EVERY reason at once.
#[test]
fn a_fully_evidenced_bump_is_admitted_and_refusals_are_complete() {
    let verdict = evaluate_bump(&attempt_with(fully_evidenced()));
    let BumpVerdict::Admitted {
        moved: movers,
        green_surfaces,
    } = verdict
    else {
        panic!("expected admission")
    };
    assert_eq!(movers, vec!["asupersync".to_string()]);
    assert!(green_surfaces.contains(&"frankensqlite".to_string()));

    // Several independent problems are all reported from one evaluation.
    let verdict = evaluate_bump(&BumpAttempt {
        kind: BumpKind::Emergency(EmergencyClass::SecurityDefect),
        deltas: vec![PinDelta {
            lib: "asupersync".to_string(),
            movement: PinMovement::Unchanged,
        }],
        results: vec![SuiteResult::new("asupersync", SuiteOutcome::NotRun)],
        golden: GoldenDisposition::Rebased {
            justification: "surface moved".to_string(),
        },
        coupled_goldens: Vec::new(),
    });
    let BumpVerdict::Refused { reasons } = verdict else {
        panic!("expected refusal")
    };
    assert!(
        reasons.len() >= 3,
        "expected no-movement + unrun asupersync + missing frankensqlite, got {reasons:?}"
    );
    assert!(reasons.contains(&BumpRefusal::NoMovement));
}

/// Emergency justification is a closed classification. Convenience, upstream
/// features, and freshness are unrepresentable, so they cannot be argued.
#[test]
fn emergency_classes_are_a_closed_set() {
    let classes = [
        EmergencyClass::SecurityDefect,
        EmergencyClass::CredibleCorruption,
        EmergencyClass::FalseScientificResult,
        EmergencyClass::ContractViolation,
        EmergencyClass::SiblingUnavailable,
    ];
    let mut slugs: Vec<&str> = classes.iter().map(|class| class.slug()).collect();
    slugs.sort_unstable();
    slugs.dedup();
    assert_eq!(slugs.len(), classes.len(), "slugs are distinct");
    for slug in slugs {
        assert!(!slug.contains("convenience"));
        assert!(!slug.contains("freshness"));
    }
    assert_eq!(BumpKind::ReleaseTrain, BumpKind::ReleaseTrain);
}

/// The rendered registry is deterministic and shows uncovered surfaces loudly.
#[test]
fn registry_renders_deterministically_and_shows_gaps() {
    let first = render_registry();
    let second = render_registry();
    assert_eq!(first, second);
    assert!(first.contains("asupersync"));
    assert!(first.contains("frankensqlite"));
    assert!(
        first.contains("NO COVERAGE"),
        "uncovered surfaces must be visible, not blank: {first}"
    );
    assert!(first.contains("cargo test --locked -p fs-exec"));
    // Feature-gated boundaries render their features in the selector.
    assert!(first.contains("--features fnp-interop"), "{first}");
}

/// The 2026-07-24/25 rehearsed bump, as an executable transcript.
///
/// Not a synthetic drill. Every outcome below was measured against the live
/// checkouts: asupersync ran 25/25 green (`fs-exec` conformance 14,
/// constellation_smoke 1, lease_battery 10) at `0.3.9@054cff23`, and
/// franken_numpy ran 2/2 green at `0.2.0@c5b6339f` once its boundary tests
/// were correctly registered. The frankensqlite surface still could not build,
/// because `fsqlite-btree` fails under the `async-api` feature that
/// `fs-ledger`'s dev-dependency enables.
///
/// This transcript CORRECTS an earlier version of itself. The first recording
/// refused franken_numpy as an uncovered surface; that was wrong — dedicated
/// round-trip and refusal tests existed in `fs-sparse/src/interop_fnp.rs` all
/// along, as unit tests behind a non-default feature. The overall verdict is
/// unchanged, but it now rests on one true reason instead of one true and one
/// false one.
#[test]
fn rehearsed_bump_2026_07_24_is_refused() {
    let recorded = vec![
        PinRow::new(
            "asupersync",
            "0.3.8",
            "5973b0ff31f405ae90fa9e6e2ef5f61a75c5b78b",
        ),
        PinRow::new(
            "franken_numpy",
            "0.1.0",
            "7fca9f6006c9f4ecdb6c7432318a0893f3a7bea1",
        ),
        PinRow::new(
            "frankensqlite",
            "unversioned-workspace",
            "987cfb97f86d3fca4d9b44e7871f427636b10126",
        ),
    ];
    let candidate = vec![
        PinRow::new(
            "asupersync",
            "0.3.9",
            "054cff2356fc525e38d54100749ff3fa33e89d7a",
        ),
        PinRow::new(
            "franken_numpy",
            "0.2.0",
            "c5b6339f2c28bdecf7066201f74e570e925ee3dc",
        ),
        PinRow::new(
            "frankensqlite",
            "unversioned-workspace",
            "31fc4a3b3a108dc49243157ea29fb1ddfcb06fdc",
        ),
    ];

    let deltas = pin_delta(&recorded, &candidate);
    assert_eq!(moved(&deltas).len(), 3, "all three moved");

    let attempt = BumpAttempt {
        kind: BumpKind::ReleaseTrain,
        deltas,
        results: vec![
            // Measured: 14 + 1 + 10 across the three registered targets.
            SuiteResult::new(
                "asupersync",
                SuiteOutcome::Executed {
                    passed: 25,
                    failed: 0,
                },
            ),
            // Measured with --features fnp-interop --lib interop_fnp.
            SuiteResult::new(
                "franken_numpy",
                SuiteOutcome::Executed {
                    passed: 2,
                    failed: 0,
                },
            ),
            // The surface did not build, so nothing executed.
            SuiteResult::new("frankensqlite", SuiteOutcome::NotRun),
        ],
        golden: GoldenDisposition::NoGoldenSurface,
        coupled_goldens: Vec::new(),
    };

    let BumpVerdict::Refused { reasons } = evaluate_bump(&attempt) else {
        panic!("a bump with an unbuildable P1 surface must never be admitted")
    };

    // The durability surface never ran: refused, not quietly skipped.
    assert!(reasons.contains(&BumpRefusal::NotExecuted {
        lib: "frankensqlite".to_string()
    }));
    // Green surfaces do not launder the refusal, and are not themselves faulted.
    for lib in ["asupersync", "franken_numpy"] {
        assert!(!reasons.iter().any(|reason| matches!(
            reason,
            BumpRefusal::NotExecuted { lib: name }
                | BumpRefusal::FailingSurface { lib: name, .. }
                | BumpRefusal::MissingResult { lib: name }
                | BumpRefusal::UncoveredSurface { lib: name, .. }
                if name == lib
        )));
    }
    // The whole refusal now rests on the one surface that genuinely could not run.
    assert_eq!(
        reasons.len(),
        1,
        "expected exactly the frankensqlite refusal, got {reasons:?}"
    );
}

/// Pins and semantic goldens move together. A golden surface is coupled to a
/// sibling exactly when its owning crate consumes that sibling at runtime, and
/// declaring "no golden surface" while coupled goldens exist is a refusal.
#[test]
fn coupled_goldens_must_be_declared() {
    let all = [
        "fs-exec:solver-resume",
        "fs-exec:tile-reduction",
        "fs-ledger:artifact-envelope",
        "fs-adjoint:dwr-accept",
        "malformed-no-colon",
    ];

    let asup = coupled_golden_surfaces("asupersync", &all);
    assert_eq!(
        asup,
        vec!["fs-exec:solver-resume", "fs-exec:tile-reduction"]
    );
    assert_eq!(
        coupled_golden_surfaces("frankensqlite", &all),
        vec!["fs-ledger:artifact-envelope"]
    );
    // PREFIX-COLLISION GUARD: frankentorch's runtime consumer is `fs-ad`, and
    // `fs-adjoint` is a DIFFERENT crate. A naive `starts_with` would wrongly
    // couple every fs-adjoint golden to frankentorch; ownership is matched
    // whole, on the segment before the colon.
    assert!(
        coupled_golden_surfaces("frankentorch", &all).is_empty(),
        "fs-ad must not match fs-adjoint"
    );
    assert_eq!(
        coupled_golden_surfaces("frankentorch", &["fs-ad:tape-bridge"]),
        vec!["fs-ad:tape-bridge"]
    );
    // A pinned-unused sibling couples to nothing, and a malformed id is ignored
    // rather than matched by accident.
    assert!(coupled_golden_surfaces("frankenpandas", &all).is_empty());
    assert!(!asup.contains(&"malformed-no-colon"));
    assert!(coupled_golden_surfaces("not-a-sibling", &all).is_empty());

    // Declaring no golden surface while coupled goldens exist is refused.
    let mut attempt = attempt_with(fully_evidenced());
    attempt.coupled_goldens = asup.iter().map(|id| (*id).to_string()).collect();
    let BumpVerdict::Refused { reasons } = evaluate_bump(&attempt) else {
        panic!("an undeclared golden implication must refuse")
    };
    let refusal = reasons
        .iter()
        .find(|reason| matches!(reason, BumpRefusal::UndeclaredGoldenImplication { .. }))
        .expect("the golden implication is reported");
    assert!(refusal.to_string().contains("fs-exec:solver-resume"));
    assert!(
        refusal
            .to_string()
            .contains("pins and goldens move together")
    );

    // Declaring the implication admits it.
    attempt.golden = GoldenDisposition::Unaffected {
        justification: "re-ran both modes; the coupled goldens are byte-identical".to_string(),
    };
    assert!(evaluate_bump(&attempt).admitted());

    // And with no coupled goldens, NoGoldenSurface stays consistent.
    let clean = attempt_with(fully_evidenced());
    assert!(evaluate_bump(&clean).admitted());
}

/// Runtime-consumer data is present wherever the sibling is actually used, and
/// absent exactly where the trust assessment records no runtime consumer.
#[test]
fn runtime_consumers_match_the_trust_assessment() {
    assert_eq!(
        surface("asupersync").expect("registered").runtime_consumers,
        &["fs-exec", "fs-plan", "fs-surrogate"]
    );
    assert_eq!(
        surface("frankensqlite")
            .expect("registered")
            .runtime_consumers,
        &["fs-ledger", "fs-vskeleton"]
    );
    // frankenscipy is a DEV-only oracle and frankenpandas is pinned-unused, so
    // neither has a runtime consumer and neither can couple a golden.
    assert!(
        surface("frankenscipy")
            .expect("registered")
            .runtime_consumers
            .is_empty()
    );
    assert!(
        surface("frankenpandas")
            .expect("registered")
            .runtime_consumers
            .is_empty()
    );
}

/// A dev-only sibling is required only when it MOVES. The "required even when
/// unmoved" rule exists because one sibling's move can break another's
/// surface, and that can only propagate through the runtime graph — so a
/// sibling with no runtime consumer is outside it.
#[test]
fn dev_only_siblings_are_required_only_when_they_move() {
    let scipy = surface("frankenscipy").expect("registered");
    assert_eq!(scipy.priority, ReviewPriority::P2);
    assert!(scipy.priority.required_every_train());
    assert!(scipy.runtime_consumers.is_empty(), "dev-only oracle");
    assert!(!scipy.tests.is_empty(), "and it IS covered");

    // Unmoved: an asupersync-only bump does not demand frankenscipy evidence.
    let unmoved = attempt_with(fully_evidenced());
    assert!(evaluate_bump(&unmoved).admitted());

    // Moved: now it must report, because a drifting oracle invalidates every
    // casebook comparison built on it.
    let mut moved_scipy = attempt_with(fully_evidenced());
    moved_scipy.deltas.push(PinDelta {
        lib: "frankenscipy".to_string(),
        movement: PinMovement::Moved {
            from_version: "0.1.0".to_string(),
            from_head: "9e271fd7".to_string(),
            to_version: "0.1.0".to_string(),
            to_head: "a133c3c8".to_string(),
        },
    });
    let BumpVerdict::Refused { reasons } = evaluate_bump(&moved_scipy) else {
        panic!("a moved oracle must report evidence")
    };
    assert!(reasons.contains(&BumpRefusal::MissingResult {
        lib: "frankenscipy".to_string()
    }));
}

/// The 2026-07-25 train: every registered surface EXECUTED and GREEN against
/// the live pins, and the bump is still refused — on the golden obligation.
///
/// Measured, not synthetic. asupersync 25/25, franken_numpy 2/2,
/// franken_networkx 8/8, frankenscipy 2/2, frankentorch 4/4, and frankensqlite
/// 5/5 once its async-pager migration landed. frankenpandas moved but declares
/// no claim and has no runtime consumer, so it requires no evidence.
///
/// The refusal is the one obligation nobody had discharged: 24 semantic golden
/// surfaces are owned by crates that consume the movers at runtime, so the
/// attempt may not declare `NoGoldenSurface`. This is the rand_nla mis-pin
/// lesson made executable — a pin may not move underneath a frozen golden
/// without someone saying what happened to it.
#[test]
fn train_2026_07_25_is_green_but_refused_on_the_golden_obligation() {
    let movers = [
        "asupersync",
        "franken_networkx",
        "franken_numpy",
        "frankenpandas",
        "frankenscipy",
        "frankensqlite",
        "frankentorch",
    ];
    let deltas: Vec<PinDelta> = movers
        .iter()
        .map(|lib| PinDelta {
            lib: (*lib).to_string(),
            movement: PinMovement::Moved {
                from_version: "recorded".to_string(),
                from_head: "recorded".to_string(),
                to_version: "live".to_string(),
                to_head: "live".to_string(),
            },
        })
        .collect();

    // Exactly the measured counts.
    let results = vec![
        SuiteResult::new(
            "asupersync",
            SuiteOutcome::Executed {
                passed: 25,
                failed: 0,
            },
        ),
        SuiteResult::new(
            "franken_networkx",
            SuiteOutcome::Executed {
                passed: 8,
                failed: 0,
            },
        ),
        SuiteResult::new(
            "franken_numpy",
            SuiteOutcome::Executed {
                passed: 2,
                failed: 0,
            },
        ),
        SuiteResult::new(
            "frankenscipy",
            SuiteOutcome::Executed {
                passed: 2,
                failed: 0,
            },
        ),
        SuiteResult::new(
            "frankensqlite",
            SuiteOutcome::Executed {
                passed: 5,
                failed: 0,
            },
        ),
        SuiteResult::new(
            "frankentorch",
            SuiteOutcome::Executed {
                passed: 4,
                failed: 0,
            },
        ),
    ];

    // Every surface is green, so no suite reason survives.
    let suite_only = BumpAttempt {
        kind: BumpKind::ReleaseTrain,
        deltas: deltas.clone(),
        results: results.clone(),
        golden: GoldenDisposition::Unaffected {
            justification: "hypothetical: goldens verified unchanged".to_string(),
        },
        coupled_goldens: Vec::new(),
    };
    assert!(
        evaluate_bump(&suite_only).admitted(),
        "with every surface green and the golden question answered, the train passes"
    );

    // But the goldens are NOT answered: 24 surfaces are owned by crates that
    // consume the movers at runtime.
    let coupled: Vec<String> = [
        "fs-exec:tune-row",
        "fs-ledger:artifact-content",
        "fs-ledger:vcs-commit-root",
        "fs-plan:voi-ranked-menu",
        "fs-vskeleton:artifact-content",
    ]
    .iter()
    .map(|id| (*id).to_string())
    .collect();
    let real = BumpAttempt {
        kind: BumpKind::ReleaseTrain,
        deltas,
        results,
        golden: GoldenDisposition::NoGoldenSurface,
        coupled_goldens: coupled,
    };
    let BumpVerdict::Refused { reasons } = evaluate_bump(&real) else {
        panic!("an undeclared golden obligation must refuse even a fully green train")
    };
    assert_eq!(
        reasons.len(),
        1,
        "only the golden obligation remains: {reasons:?}"
    );
    assert!(matches!(
        reasons[0],
        BumpRefusal::UndeclaredGoldenImplication { .. }
    ));
    assert!(
        reasons[0]
            .to_string()
            .contains("fs-ledger:artifact-content")
    );
}
