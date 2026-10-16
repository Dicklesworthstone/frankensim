//! Binding-witness algebra battery (bead frankensim-h0vur, core slice):
//! the meet's argmin must survive composition — typed, nonempty, canonical,
//! duplicate-free, tie-complete, and impossible to launder. The spine e2e
//! discrimination test named by the bead stays with the bead (non-vacuous
//! only once conduction/QoI execute); this battery proves the algebra.

use fs_evidence::ValidityDomain;
use fs_evidence::color::{Color, IntervalOp};
use fs_evidence::witness::{
    BindingCause, WitnessError, WitnessedColor, compose_all_witnessed, compose_witnessed,
};

fn verified() -> Color {
    Color::Verified { lo: 1.0, hi: 2.0 }
}

fn validated(regime: ValidityDomain, dataset: &str) -> Color {
    Color::Validated {
        regime,
        dataset: dataset.to_string(),
    }
}

fn estimated(estimator: &str) -> Color {
    Color::Estimated {
        estimator: estimator.to_string(),
        dispersion: 0.1,
    }
}

fn leaf(color: Color, stage: &str) -> WitnessedColor {
    WitnessedColor::leaf(color, stage, "").expect("nonempty stage")
}

fn ids(w: &WitnessedColor) -> Vec<String> {
    w.witnesses().iter().map(|b| b.id.to_string()).collect()
}

#[test]
fn every_rank_pair_attributes_the_weakest_operand() {
    let re = ValidityDomain::unconstrained().with("Re", 1e3, 1e5);
    // (a colour, its rank label) fixtures.
    let cases: Vec<(Color, &str)> = vec![
        (verified(), "verified-stage"),
        (validated(re.clone(), "wind-tunnel"), "validated-stage"),
        (estimated("surrogate"), "estimated-stage"),
    ];
    for (ca, sa) in &cases {
        for (cb, sb) in &cases {
            let a = leaf(ca.clone(), sa);
            let b = leaf(cb.clone(), sb);
            let composed = compose_witnessed(&a, &b, IntervalOp::Add).expect("composes");
            let rank = composed.color().rank();
            // Anti-laundering: every witness's originating leaf rank must
            // equal the composed rank (no degradation is possible in these
            // compatible-regime pairs).
            for witness in composed.witnesses() {
                let origin_rank = cases
                    .iter()
                    .find(|(_, s)| witness.id.stage == **s)
                    .map(|(c, _)| c.rank())
                    .expect("witness maps to a fixture stage");
                assert_eq!(
                    origin_rank, rank,
                    "witness {witness:?} outranks composed {rank:?}"
                );
                assert_eq!(witness.cause, BindingCause::WeakestOperand);
            }
            // Ties retain BOTH; strict orders retain exactly the weaker.
            let expected = usize::from(ca.rank() == cb.rank() && sa != sb) + 1;
            assert_eq!(
                composed.witnesses().len(),
                expected,
                "{sa} + {sb}: witness cardinality"
            );
        }
    }
}

#[test]
fn discrimination_same_colour_different_stage_different_witness() {
    // The test that would have caught the ProvenanceClass defect: identical
    // composed colour, different binding stage, MUST differ in witness.
    let strong = leaf(verified(), "assign");
    let weak_a = leaf(estimated("correlation-card-x"), "flow-network");
    let weak_b = leaf(estimated("correlation-card-x"), "conduction");
    let run_a = compose_witnessed(&strong, &weak_a, IntervalOp::Add).expect("composes");
    let run_b = compose_witnessed(&strong, &weak_b, IntervalOp::Add).expect("composes");
    assert_eq!(run_a.color(), run_b.color(), "identical composed colour");
    assert_ne!(
        ids(&run_a),
        ids(&run_b),
        "different binding stage must show"
    );
    assert_eq!(ids(&run_a), vec!["flow-network".to_string()]);
    assert_eq!(ids(&run_b), vec!["conduction".to_string()]);
}

#[test]
fn multiway_ties_are_permutation_and_grouping_invariant() {
    let stages = ["s1", "s2", "s3", "s4"];
    let leaves: Vec<WitnessedColor> = stages.iter().map(|s| leaf(estimated("est"), s)).collect();
    let forward = compose_all_witnessed(&leaves, IntervalOp::Add).expect("composes");
    let mut reversed_items = leaves.clone();
    reversed_items.reverse();
    let reversed = compose_all_witnessed(&reversed_items, IntervalOp::Add).expect("composes");
    // Canonical order makes permutation invisible.
    assert_eq!(ids(&forward), ids(&reversed));
    assert_eq!(ids(&forward), vec!["s1", "s2", "s3", "s4"]);
    // Grouping invariance: ((s1+s2)+(s3+s4)) == fold order.
    let left = compose_witnessed(&leaves[0], &leaves[1], IntervalOp::Add).expect("composes");
    let right = compose_witnessed(&leaves[2], &leaves[3], IntervalOp::Add).expect("composes");
    let grouped = compose_witnessed(&left, &right, IntervalOp::Add).expect("composes");
    assert_eq!(ids(&grouped), ids(&forward));
}

#[test]
fn duplicate_witnesses_collapse() {
    let shared = leaf(estimated("est"), "flow-network");
    let composed = compose_witnessed(&shared, &shared.clone(), IntervalOp::Add).expect("composes");
    assert_eq!(ids(&composed), vec!["flow-network".to_string()]);
}

#[test]
fn disjoint_validated_regimes_degrade_with_both_witnesses() {
    // compose() lands BELOW both operands here; blaming one side would be
    // dishonest, so both witnesses survive with the degradation cause.
    let low = ValidityDomain::unconstrained().with("Re", 1.0, 10.0);
    let high = ValidityDomain::unconstrained().with("Re", 100.0, 1000.0);
    let a = leaf(validated(low, "tunnel-low"), "airflow");
    let b = leaf(validated(high, "tunnel-high"), "convection");
    let composed = compose_witnessed(&a, &b, IntervalOp::Add).expect("composes");
    assert_eq!(
        composed.color().rank(),
        fs_evidence::ColorRank::Estimated,
        "disjoint regimes demote"
    );
    assert_eq!(ids(&composed), vec!["airflow", "convection"]);
    assert!(
        composed
            .witnesses()
            .iter()
            .all(|w| w.cause == BindingCause::CompositionDegraded),
        "degradation is the composition's doing, both sides carry the cause"
    );
}

#[test]
fn cardinality_and_empty_stage_refuse_by_name() {
    let err = WitnessedColor::leaf(estimated("e"), "", "").unwrap_err();
    assert!(err.to_string().contains("FS-EVIDENCE-WITNESS-EMPTY-STAGE"));

    let leaves: Vec<WitnessedColor> = (0..70)
        .map(|i| leaf(estimated("e"), &format!("stage-{i:03}")))
        .collect();
    let err = compose_all_witnessed(&leaves, IntervalOp::Add).unwrap_err();
    assert!(
        matches!(err, WitnessError::Cardinality { attempted } if attempted > 64),
        "{err}"
    );
    assert!(err.to_string().contains("FS-EVIDENCE-WITNESS-CARDINALITY"));
}

#[test]
fn detail_discriminant_changes_witness_at_constant_stage_and_colour() {
    let strong = leaf(verified(), "assign");
    let card_a = WitnessedColor::leaf(estimated("est"), "flow-network", "card-a").expect("leaf");
    let card_b = WitnessedColor::leaf(estimated("est"), "flow-network", "card-b").expect("leaf");
    let run_a = compose_witnessed(&strong, &card_a, IntervalOp::Add).expect("composes");
    let run_b = compose_witnessed(&strong, &card_b, IntervalOp::Add).expect("composes");
    assert_eq!(run_a.color(), run_b.color());
    assert_eq!(ids(&run_a), vec!["flow-network:card-a"]);
    assert_eq!(ids(&run_b), vec!["flow-network:card-b"]);
}
