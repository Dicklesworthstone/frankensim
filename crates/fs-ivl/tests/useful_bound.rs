//! G0/G3 driver integration for contextual useful bounds.

use fs_evidence::{ClaimClass, NoUsefulBoundCause, UsefulnessCriterion};
use fs_ivl::{Interval, RootSearchConfig, TaylorModel1, newton_roots_bounded};

fn criterion(max_width: f64, unit: &str) -> UsefulnessCriterion {
    UsefulnessCriterion::try_new("finite-horizon thermal trip", unit, max_width)
        .expect("valid criterion")
}

#[test]
fn taylor_driver_classifies_finite_and_overflowed_enclosures() {
    let domain = Interval::new(-1.0, 1.0);
    let constant = TaylorModel1::constant(4.0, domain, 1).expect("constant model");
    let useful = constant
        .bound_with_usefulness(
            criterion(1.0, "kelvin"),
            NoUsefulBoundCause::HorizonTooLong,
            ClaimClass::LongHorizonMeanLoad,
        )
        .expect("useful finite bound");
    assert!(useful.bound().is_some());

    let overflowed = TaylorModel1::constant(f64::MAX, domain, 1)
        .expect("finite constant")
        .scale(2.0)
        .expect("overflow degrades to a whole enclosure");
    let refused = overflowed
        .bound_with_usefulness(
            criterion(1.0, "kelvin"),
            NoUsefulBoundCause::HorizonTooLong,
            ClaimClass::LongHorizonMeanLoad,
        )
        .expect("whole enclosure is a typed outcome");
    assert_eq!(
        refused
            .no_useful_bound()
            .expect("whole enclosure is not useful")
            .cause(),
        NoUsefulBoundCause::LipschitzBlowup
    );
}

#[test]
fn budget_exhausted_root_driver_logs_width_trajectory_and_refuses() {
    let domain = Interval::new(-1.0, 1.0);
    let zero = |_x: Interval| Interval::point(0.0);
    let report = newton_roots_bounded(
        &zero,
        &zero,
        domain,
        RootSearchConfig {
            min_width: f64::MIN_POSITIVE,
            max_boxes: 3,
        },
    )
    .expect("bounded search");

    assert!(!report.complete);
    assert_eq!(report.width_trajectory.len(), report.boxes_examined);
    assert_eq!(report.width_trajectory[0], domain.width());
    assert!(
        report
            .width_trajectory
            .windows(2)
            .all(|pair| pair[1] <= pair[0])
    );

    let outcome = report
        .bound_with_usefulness(
            criterion(100.0, "seconds"),
            NoUsefulBoundCause::HorizonTooLong,
            ClaimClass::RootOrEventTime,
        )
        .expect("projection")
        .expect("ambiguous roots retain a hull");
    let refusal = outcome.no_useful_bound().expect("budget is absorbing");
    assert_eq!(refusal.cause(), NoUsefulBoundCause::BudgetExhausted);
    assert_eq!(
        refusal.suggested_reformulation(),
        ClaimClass::RootOrEventTime
    );
}
