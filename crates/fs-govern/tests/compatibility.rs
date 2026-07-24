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
    PinMovement, PinRow, ReviewPriority, SURFACES, SuiteOutcome, SuiteResult, evaluate_bump, moved,
    pin_delta, render_registry, surface,
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
                entry.selector().is_none(),
                "{} must not render a runnable selector with no tests",
                entry.lib
            );
        } else {
            assert!(
                entry.no_test_reason.is_none(),
                "{} has tests, so a no-test reason is contradictory",
                entry.lib
            );
            let selector = entry.selector().expect("covered surface has a selector");
            assert!(selector.starts_with("cargo test --locked"));
            assert!(selector.contains("-p "));
            assert!(selector.contains("--test "));
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

/// A moved sibling with no compatibility coverage is refused rather than
/// waved through: an uncovered surface cannot supply evidence.
#[test]
fn moved_but_uncovered_sibling_is_refused() {
    let mut attempt = attempt_with(vec![green("asupersync"), green("frankensqlite")]);
    attempt.deltas.push(PinDelta {
        lib: "frankentorch".to_string(),
        movement: PinMovement::Moved {
            from_version: "0.1.0".to_string(),
            from_head: "f00c3ce".to_string(),
            to_version: "0.2.0".to_string(),
            to_head: "deadbee".to_string(),
        },
    });
    let BumpVerdict::Refused { reasons } = evaluate_bump(&attempt) else {
        panic!("expected refusal")
    };
    let uncovered = reasons
        .iter()
        .find(|reason| matches!(reason, BumpRefusal::UncoveredSurface { lib, .. } if lib == "frankentorch"))
        .expect("uncovered surface is reported");
    assert!(uncovered.to_string().contains("feature-gated"));

    // A sibling nobody registered cannot be adjudicated at all.
    let mut unknown = attempt_with(vec![green("asupersync"), green("frankensqlite")]);
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
        results: vec![green("asupersync"), green("frankensqlite")],
        golden: GoldenDisposition::NoGoldenSurface,
    };
    let BumpVerdict::Refused { reasons } = evaluate_bump(&attempt) else {
        panic!("expected refusal")
    };
    assert!(reasons.contains(&BumpRefusal::NoMovement));
}

/// The admitted path, and the fact that refusal reports EVERY reason at once.
#[test]
fn a_fully_evidenced_bump_is_admitted_and_refusals_are_complete() {
    let verdict = evaluate_bump(&attempt_with(vec![
        green("asupersync"),
        green("frankensqlite"),
    ]));
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
}

/// The 2026-07-24 rehearsed bump, as an executable transcript.
///
/// This is not a synthetic drill. It records the real seven-sibling drift,
/// the real executed outcomes — the asupersync surface ran 25/25 green
/// (conformance 14, constellation_smoke 1, lease_battery 10) against
/// `0.3.9@054cff23`, while the frankensqlite surface could not even build
/// because `fsqlite-btree` is mid async-pager migration — and the verdict the
/// protocol actually returns. A surface that fails to compile is `NotRun`:
/// there is no execution to report, and an unbuildable dependency must never
/// read as an absent problem.
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
        // A MINOR version move, into a surface with no coverage at all.
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
            // The surface did not build, so nothing executed.
            SuiteResult::new("frankensqlite", SuiteOutcome::NotRun),
        ],
        golden: GoldenDisposition::NoGoldenSurface,
    };

    let BumpVerdict::Refused { reasons } = evaluate_bump(&attempt) else {
        panic!("a bump with an unbuildable P1 surface must never be admitted")
    };

    // The durability surface never ran: refused, not quietly skipped.
    assert!(reasons.contains(&BumpRefusal::NotExecuted {
        lib: "frankensqlite".to_string()
    }));
    // franken_numpy moved a minor version into a surface with no coverage.
    assert!(reasons.iter().any(|reason| matches!(
        reason,
        BumpRefusal::UncoveredSurface { lib, .. } if lib == "franken_numpy"
    )));
    // A green surface does not launder the refusal for the others.
    assert!(!reasons.iter().any(|reason| matches!(
        reason,
        BumpRefusal::NotExecuted { lib } | BumpRefusal::FailingSurface { lib, .. }
            if lib == "asupersync"
    )));
}
