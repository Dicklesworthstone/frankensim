//! Horizon trigger 4 battery (bead frankensim-epic-addendum-xpck.5.4):
//! boundary, mutant, and property coverage for the splitting-error demand
//! gate in `fs_govern::horizon_splitting`.

use fs_govern::horizon_splitting::{
    current_receipt, population_disposition, splitting_verdict, BudgetRefusal, ErrorTerm,
    PopulationDisposition, SplittingVerdict, WorkloadBudget, WorkloadClass, SPLITTING_SHARE_MIN,
};

fn budget(class: WorkloadClass, terms: &[(&str, f64)], stable: bool) -> WorkloadBudget {
    WorkloadBudget {
        workload_id: "wl-1".into(),
        class,
        terms: terms.iter().map(|(n, s)| ErrorTerm::new(n, *s)).collect(),
        coupling_stable: stable,
    }
}

fn paying_dominant() -> WorkloadBudget {
    budget(
        WorkloadClass::Paying,
        &[("splitting", 0.45), ("discretization", 0.25), ("iteration", 0.2), ("model", 0.1)],
        true,
    )
}

#[test]
fn gate_fires_only_for_paying_stable_strictly_dominant_at_or_above_threshold() {
    assert_eq!(splitting_verdict(&paying_dominant()), Ok(SplittingVerdict::Activate));

    // Same shape but non-paying: measured, never activated.
    let mut w = paying_dominant();
    w.class = WorkloadClass::NonPaying;
    assert_eq!(splitting_verdict(&w), Ok(SplittingVerdict::InstrumentOnly));

    // Unstable coupling: premise gone, nonactivating even when dominant.
    let mut w = paying_dominant();
    w.coupling_stable = false;
    assert_eq!(splitting_verdict(&w), Ok(SplittingVerdict::InstrumentOnly));
}

#[test]
fn boundary_equality_with_threshold_activates_but_tie_does_not() {
    // Exactly at threshold AND strictly largest: fires.
    let at = budget(
        WorkloadClass::Paying,
        &[("splitting", SPLITTING_SHARE_MIN), ("discretization", 0.3), ("iteration", 0.3), ("model", 0.15)],
        true,
    );
    assert_eq!(splitting_verdict(&at), Ok(SplittingVerdict::Activate));

    // Strict tie with the largest competitor: NOT dominant, does not fire
    // (equality is not dominance).
    let tied = budget(
        WorkloadClass::Paying,
        &[("splitting", 0.35), ("discretization", 0.35), ("iteration", 0.2), ("model", 0.1)],
        true,
    );
    assert_eq!(splitting_verdict(&tied), Ok(SplittingVerdict::InstrumentOnly));
}

#[test]
fn malformed_budgets_refuse_by_name_never_verdict() {
    // Zero budget.
    let empty = budget(WorkloadClass::Paying, &[], true);
    assert_eq!(splitting_verdict(&empty), Err(BudgetRefusal::EmptyTerms));

    // Missing required term.
    let missing = budget(
        WorkloadClass::Paying,
        &[("splitting", 0.5), ("discretization", 0.5)],
        true,
    );
    assert_eq!(
        splitting_verdict(&missing),
        Err(BudgetRefusal::MissingRequiredTerm { term: "iteration".into() })
    );

    // Double-counted error (shares sum > 1) refuses.
    let doubled = budget(
        WorkloadClass::Paying,
        &[("splitting", 0.6), ("discretization", 0.4), ("iteration", 0.2), ("model", 0.1)],
        true,
    );
    assert!(matches!(splitting_verdict(&doubled), Err(BudgetRefusal::SharesDontSum { .. })));

    // Out-of-range share refuses.
    let negative = budget(
        WorkloadClass::Paying,
        &[("splitting", -0.1), ("discretization", 0.6), ("iteration", 0.3), ("model", 0.2)],
        true,
    );
    assert!(matches!(
        splitting_verdict(&negative),
        Err(BudgetRefusal::ShareOutOfRange { .. })
    ));
}

#[test]
fn mutant_drop_larger_competitor_must_flip_the_verdict_not_hide_it() {
    // With discretization present at 0.45 vs splitting 0.40: no fire.
    let honest = budget(
        WorkloadClass::Paying,
        &[("splitting", 0.40), ("discretization", 0.45), ("iteration", 0.1), ("model", 0.05)],
        true,
    );
    assert_eq!(splitting_verdict(&honest), Ok(SplittingVerdict::InstrumentOnly));

    // MUTANT: dropping the larger competing term would make splitting
    // dominant. The engine must produce a DIFFERENT verdict for that
    // different budget — i.e. the ranking provably reads every term.
    let mutated = budget(
        WorkloadClass::Paying,
        &[("splitting", 0.40), ("iteration", 0.35), ("model", 0.25)],
        true,
    );
    assert_eq!(
        splitting_verdict(&mutated).expect("mutated budget is complete"),
        SplittingVerdict::Activate,
        "verdict must track the actual term set (no silent term dropping)"
    );
    // And the honest budget stays InstrumentOnly — the pair proves the
    // comparator consumed the dropped term rather than ignoring it.
    assert_eq!(splitting_verdict(&honest), Ok(SplittingVerdict::InstrumentOnly));
}

#[test]
fn property_renormalizing_a_uniform_scaling_preserves_the_verdict() {
    let base = paying_dominant();
    let expected = splitting_verdict(&base);
    for k in [0.5_f64, 2.0, 10.0] {
        let scaled_terms: Vec<ErrorTerm> = base
            .terms
            .iter()
            .map(|t| ErrorTerm::new(&t.name, t.share * k))
            .collect();
        let sum: f64 = scaled_terms.iter().map(|t| t.share).sum();
        let renorm: Vec<ErrorTerm> = scaled_terms
            .into_iter()
            .map(|t| ErrorTerm::new(&t.name, t.share / sum))
            .collect();
        let scaled = WorkloadBudget { terms: renorm, ..base.clone() };
        assert_eq!(
            splitting_verdict(&scaled),
            expected,
            "unit rescale + renormalization must preserve ratios and verdicts"
        );
    }
}

#[test]
fn property_workload_order_never_changes_the_population_class() {
    let a = paying_dominant();
    let b = budget(
        WorkloadClass::Paying,
        &[("splitting", 0.15), ("discretization", 0.45), ("iteration", 0.25), ("model", 0.15)],
        true,
    );
    let c = budget(
        WorkloadClass::NonPaying,
        &[("splitting", 0.5), ("discretization", 0.2), ("iteration", 0.2), ("model", 0.1)],
        true,
    );
    let forward = population_disposition(&[a.clone(), b.clone(), c.clone()]);
    let backward = population_disposition(&[c, b, a]);
    assert_eq!(
        forward, backward,
        "population disposition must be order-independent"
    );
    assert!(matches!(forward, PopulationDisposition::InstrumentOnly { .. }));
}

#[test]
fn empty_and_nonpaying_populations_yield_nodata() {
    assert!(matches!(
        population_disposition(&[]),
        PopulationDisposition::NoData { .. }
    ));
    let only_nonpaying =
        budget(WorkloadClass::NonPaying, &[("splitting", 0.9), ("discretization", 0.05), ("iteration", 0.03), ("model", 0.02)], true);
    assert!(matches!(
        population_disposition(&[only_nonpaying]),
        PopulationDisposition::NoData { .. }
    ));
}

#[test]
fn all_paying_passing_population_activates_weakest_link_holds() {
    let a = paying_dominant();
    let b = budget(
        WorkloadClass::Paying,
        &[("splitting", 0.30), ("discretization", 0.28), ("iteration", 0.24), ("model", 0.18)],
        true,
    );
    assert_eq!(
        population_disposition(&[a.clone(), b]),
        PopulationDisposition::Activate {
            verdicts: vec![
                ("wl-1".to_string(), SplittingVerdict::Activate),
                ("wl-1".to_string(), SplittingVerdict::Activate),
            ],
        }
    );
    // One refusing budget blocks activation for everyone.
    let mut broken = paying_dominant();
    broken.workload_id = "broken".into();
    broken.terms.clear();
    match population_disposition(&[a, broken]) {
        PopulationDisposition::InstrumentOnly { verdicts } => {
            assert!(verdicts.iter().any(|(id, _)| id.contains("refused")));
        }
        other => panic!("expected InstrumentOnly, got {other:?}"),
    }
}

#[test]
fn today_is_nodata_and_that_is_typed_not_assumed() {
    match current_receipt() {
        PopulationDisposition::NoData { reason } => {
            assert!(!reason.is_empty());
        }
        other => panic!("program has no paying coupled-transient workload yet: {other:?}"),
    }
}
