//! G0 policy tests for the fixed-size moonshot portfolio.

use fs_govern::{
    DisplacementRecord, MOONSHOT_V1_INITIAL_CAP, MoonshotBudget, MoonshotDeclaration,
    MoonshotDisposition, MoonshotError, MoonshotPortfolio, NamedFalsifier, QuarterlyReview,
    ReplacementAdmission,
};

fn declaration(
    bead_id: &str,
    lane: &str,
    legacy_baseline: bool,
    effort_minutes: u64,
) -> MoonshotDeclaration {
    MoonshotDeclaration::new(
        bead_id,
        lane,
        "test-owner",
        legacy_baseline,
        NamedFalsifier::new(
            "counterexample",
            "the registered counterexample violates the proposed mechanism",
            "terminalize the mechanism and retain the counterexample",
        )
        .expect("falsifier"),
        MoonshotBudget::new(effort_minutes, "2026-10-01").expect("budget"),
        "the moonshot receives no product-critical-path capacity",
        QuarterlyReview::new(
            "2026-10-01",
            "scheduled",
            "independently checked counterexample corpus",
        )
        .expect("review"),
    )
    .expect("declaration")
}

fn displacement(bead_id: &str, disposition: MoonshotDisposition) -> DisplacementRecord {
    DisplacementRecord::new(
        bead_id,
        disposition,
        "artifact:test-state",
        "bounded test disposition",
    )
    .expect("displacement")
}

#[test]
fn moonshot_declarations_refuse_unbounded_or_incoherent_policy_fields() {
    assert_eq!(
        MoonshotBudget::new(0, "2026-10-01"),
        Err(MoonshotError::ZeroEffortBudget)
    );
    assert_eq!(
        MoonshotBudget::new(1, "2026-02-30"),
        Err(MoonshotError::InvalidDate {
            field: "budget.calendar_deadline"
        })
    );

    let late_review = MoonshotDeclaration::new(
        "bead-a",
        "lane-a",
        "owner",
        true,
        NamedFalsifier::new("f", "observation", "decision").expect("falsifier"),
        MoonshotBudget::new(60, "2026-09-30").expect("budget"),
        "off path",
        QuarterlyReview::new("2026-10-01", "scheduled", "evidence").expect("review"),
    );
    assert_eq!(late_review, Err(MoonshotError::ReviewAfterDeadline));

    let over_cap = MoonshotPortfolio::new(MOONSHOT_V1_INITIAL_CAP + 1, Vec::new());
    assert!(matches!(over_cap, Err(MoonshotError::CapExceeded { .. })));
}

#[test]
fn moonshot_one_for_one_replacement_releases_exactly_one_slot_and_retains_disposition() {
    let portfolio = MoonshotPortfolio::new(
        2,
        vec![
            declaration("bead-a", "lane-a", true, 100),
            declaration("bead-b", "lane-b", true, 100),
        ],
    )
    .expect("baseline portfolio");
    let replacement = ReplacementAdmission::new(
        declaration("bead-c", "lane-c", false, 80),
        displacement("bead-a", MoonshotDisposition::ShelvedWithState),
    )
    .expect("replacement request");
    let portfolio = portfolio
        .assess_replacement(&replacement)
        .expect("admissible replacement");

    assert_eq!(portfolio.cap(), 2);
    assert_eq!(portfolio.active_count(), 2);
    assert!(!portfolio.is_active("bead-a"));
    assert!(portfolio.is_terminal("bead-a"));
    assert!(portfolio.is_active("bead-b"));
    assert!(portfolio.is_active("bead-c"));
    assert_eq!(portfolio.displacements().len(), 1);
    assert_eq!(
        portfolio
            .displacements()
            .first()
            .expect("one displacement")
            .disposition(),
        MoonshotDisposition::ShelvedWithState
    );
    assert_eq!(portfolio.ledger().decisions().len(), 4);

    let revive = ReplacementAdmission::new(
        declaration("bead-a", "lane-a-revival", false, 20),
        displacement("bead-c", MoonshotDisposition::Completed),
    )
    .expect("syntactically complete revival request");
    assert!(matches!(
        portfolio.assess_replacement(&revive),
        Err(MoonshotError::TerminalBeadRevival { bead_id }) if bead_id == "bead-a"
    ));
}

#[test]
fn moonshot_replacement_refuses_legacy_candidates_and_budget_growth() {
    let legacy_candidate = ReplacementAdmission::new(
        declaration("bead-c", "lane-c", true, 50),
        displacement("bead-a", MoonshotDisposition::Completed),
    );
    assert!(matches!(
        legacy_candidate,
        Err(MoonshotError::CandidateIsLegacyBaseline)
    ));

    let portfolio = MoonshotPortfolio::new(
        2,
        vec![
            declaration("bead-a", "lane-a", true, 100),
            declaration("bead-b", "lane-b", true, 100),
        ],
    )
    .expect("baseline portfolio");
    let over_budget = ReplacementAdmission::new(
        declaration("bead-c", "lane-c", false, 201),
        displacement("bead-a", MoonshotDisposition::Falsified),
    )
    .expect("replacement request");
    assert!(matches!(
        portfolio.assess_replacement(&over_budget),
        Err(MoonshotError::Lane(_))
    ));
}

#[test]
fn moonshot_liquidation_shrinks_the_cap_without_leaving_a_refillable_slot() {
    let portfolio = MoonshotPortfolio::new(
        2,
        vec![
            declaration("bead-a", "lane-a", true, 100),
            declaration("bead-b", "lane-b", true, 100),
        ],
    )
    .expect("baseline portfolio")
    .liquidate(&displacement("bead-a", MoonshotDisposition::Completed))
    .expect("liquidation");
    assert_eq!(portfolio.cap(), 1);
    assert_eq!(portfolio.active_count(), 1);
    assert!(portfolio.is_terminal("bead-a"));

    let replacement = ReplacementAdmission::new(
        declaration("bead-c", "lane-c", false, 100),
        displacement("bead-b", MoonshotDisposition::Falsified),
    )
    .expect("replacement request");
    let portfolio = portfolio
        .assess_replacement(&replacement)
        .expect("replacement at shrunken cap");
    assert_eq!(portfolio.cap(), 1);
    assert_eq!(portfolio.active_count(), 1);
    assert!(portfolio.is_active("bead-c"));
}
