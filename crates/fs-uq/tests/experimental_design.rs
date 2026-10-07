use fs_uq::experimental_design::*;
fn group(name: &str, rows: &[[f64; 2]]) -> ExperimentGroup {
    ExperimentGroup { name: name.into(), rows: rows.iter().map(|r| r.to_vec()).collect() }
}
fn config(k: usize) -> ExperimentDesignConfig {
    ExperimentDesignConfig { parameters: 2, maximum_selected: k, ridge_precision: 1.0,
        maximum_factorizations: 4096 }
}
fn candidates() -> Vec<ExperimentGroup> {
    vec![group("strong-x", &[[10.0, 0.0]]), group("redundant-x", &[[9.0, 0.0]]),
        group("independent-y", &[[0.0, 4.0]])]
}
#[test]
fn complementary_experiments_beat_redundant_high_sensitivity() {
    let groups = candidates(); let before = groups.clone();
    let r = select_experiments(&groups, config(2), || false).unwrap();
    assert_eq!(r.selected, vec![0, 2]);
    assert_eq!(r.stop, ExperimentDesignStop::SelectionLimit);
    assert_eq!(r.factorizations, 6);
    assert!((r.log_determinant_gain - (101.0_f64*17.0).ln()).abs() < 1e-12);
    assert!((r.steps[1].marginal_gain - 17.0_f64.ln()).abs() < 1e-12);
    assert_eq!(groups, before);
    assert_eq!(r, select_experiments(&groups, config(2), || false).unwrap());
}
#[test]
fn whole_blocks_are_indivisible_and_ties_follow_declaration_order() {
    let groups = vec![group("one", &[[3.0, 0.0], [0.0, 3.0]]),
        group("two", &[[3.0, 0.0], [0.0, 3.0]]), group("one-large-row", &[[9.0, 0.0]])];
    let r = select_experiments(&groups, config(1), || false).unwrap();
    assert_eq!(r.selected, vec![0]); // det=100, not the row's det=82
    assert!((r.log_determinant_gain - 100.0_f64.ln()).abs() < 1e-12);
}
#[test]
fn incomplete_score_round_never_installs_an_input_order_winner() {
    let groups = candidates();
    let mut c = config(2); c.maximum_factorizations = 3;
    let r = select_experiments(&groups, c, || false).unwrap();
    assert_eq!(r.stop, ExperimentDesignStop::EvaluationLimit);
    assert!(r.selected.is_empty()); assert_eq!(r.factorizations, 1);
    c.maximum_factorizations = 5;
    let r = select_experiments(&groups, c, || false).unwrap();
    assert_eq!(r.stop, ExperimentDesignStop::EvaluationLimit);
    assert_eq!(r.selected, vec![0]); assert_eq!(r.factorizations, 4);
}
#[test]
fn zero_rows_do_not_claim_information_or_force_an_arbitrary_selection() {
    let groups = vec![group("zero", &[[0.0, 0.0]]), group("signal", &[[1.0, 0.0]])];
    let r = select_experiments(&groups, config(2), || false).unwrap();
    assert_eq!(r.selected, vec![1]); assert_eq!(r.stop, ExperimentDesignStop::NoResolvedGain);
    assert!((r.log_determinant_gain - 2.0_f64.ln()).abs() < 1e-12);
}
#[test]
fn common_jacobian_and_precision_rescaling_preserve_the_design() {
    let expected = select_experiments(&candidates(), config(2), || false).unwrap();
    for scale in [1e-100, 0.01, 100.0, 1e100] {
        let mut groups = candidates(); let mut c = config(2); c.ridge_precision = scale*scale;
        for g in &mut groups { for r in &mut g.rows { for x in r { *x *= scale; } } }
        let actual = select_experiments(&groups, c, || false).unwrap();
        assert_eq!(actual.selected, expected.selected);
        assert!((actual.log_determinant_gain - expected.log_determinant_gain).abs() < 1e-12);
    }
}
#[test]
fn invalid_families_and_numerical_ranges_are_refused() {
    for value in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        let mut c = config(1); c.ridge_precision = value;
        assert!(select_experiments(&candidates(), c, || false).is_err());
    }
    assert!(select_experiments(&[], config(1), || false).is_err());
    let mut groups = candidates(); groups[1].name = groups[0].name.clone();
    assert!(select_experiments(&groups, config(1), || false).is_err());
    for rows in [vec![], vec![vec![1.0]], vec![vec![f64::NAN, 0.0]], vec![vec![f64::MAX, 0.0]]] {
        let mut groups = candidates(); groups[2].rows = rows;
        assert!(select_experiments(&groups, config(1), || false).is_err());
    }
    let mut c = config(1); c.parameters = usize::MAX;
    assert!(select_experiments(&candidates(), c, || false).is_err());
    c = config(4); assert!(select_experiments(&candidates(), c, || false).is_err());
}
#[test]
fn cancellation_at_every_boundary_returns_no_partial_selection() {
    let groups = candidates(); let mut polls = 0;
    let expected = select_experiments(&groups, config(2), || { polls += 1; false }).unwrap();
    for stop in 1..=polls {
        let mut seen = 0;
        assert_eq!(select_experiments(&groups, config(2), || { seen += 1; seen == stop }),
            Err(ExperimentDesignError::Cancelled));
    }
    assert_eq!(select_experiments(&groups, config(2), || false).unwrap(), expected);
}
