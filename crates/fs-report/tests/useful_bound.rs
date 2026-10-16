//! Honest visual rendering for `NoUsefulBound`.

use fs_evidence::{
    BoundInterval, BoundOutcome, ClaimClass, NoUsefulBoundCause, UsefulnessCriterion,
};
use fs_report::no_useful_bound_markdown;

#[test]
fn report_uses_a_distinct_visual_class_and_actionable_reformulation() {
    let criterion =
        UsefulnessCriterion::try_new("chaotic thermal trajectory at 24 hours", "kelvin", 0.5)
            .expect("criterion");
    let outcome = BoundOutcome::refuse(
        BoundInterval::try_new(f64::NEG_INFINITY, f64::INFINITY).expect("whole enclosure"),
        criterion,
        NoUsefulBoundCause::LipschitzBlowup,
        ClaimClass::LongHorizonMeanLoad,
    );
    let markdown =
        no_useful_bound_markdown(outcome.no_useful_bound().expect("typed usefulness refusal"));

    assert!(markdown.starts_with("### NoUsefulBound\n"));
    assert!(markdown.contains("`lipschitz-blowup`"));
    assert!(markdown.contains("`long-horizon-mean-load`"));
    assert!(markdown.contains("Long-horizon mean load"));
    assert!(markdown.contains("no compliance verdict, scientific color, or certificate"));
}
