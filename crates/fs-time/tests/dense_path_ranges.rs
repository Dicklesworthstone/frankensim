use fs_time::adaptive::events::multiple::{
    EventOrder, EventSetAdvance, EventSetOptions, EventSpec, rk45_until_any_event,
    rk45_until_any_event_observed,
};
use fs_time::{
    AdaptiveState, AdaptiveStatus, DenseAdvance, EventDirection, EventError, InitialEvent,
    PiController, rk45_dense_step,
};
use std::cell::Cell;

fn options() -> EventSetOptions {
    EventSetOptions {
        max_step: 1.0,
        scan_substeps: 8,
        time_tolerance: 1e-12,
        max_iterations: 80,
        order: EventOrder::RequireSeparated,
    }
}

#[test]
fn range_includes_internal_extremum_and_restricts_the_accepted_prefix() {
    let rhs = |t: f64, _: &[f64], out: &mut [f64]| {
        out[0] = 4.0 - 8.0 * t;
        out[1] = 0.0;
    };
    let mut state = AdaptiveState::new(0.0, &[0.0, -3.0], 1.0);
    let DenseAdvance::Accepted(dense) = rk45_dense_step(
        &mut state,
        &rhs,
        1.0,
        1e-11,
        1e-12,
        &PiController::default(),
        4,
        &mut || false,
    )
    .unwrap() else {
        panic!("the exact quadratic should be accepted");
    };
    let full = dense.component_range(0, 1.0).unwrap();
    assert!(full[0] <= 0.0 && full[1] >= 1.0 && full[1] < 1.34);
    for i in 0..=1000 {
        let value = dense.evaluate(i as f64 / 1000.0).unwrap()[0];
        assert!(value >= full[0] && value <= full[1]);
    }
    let prefix = dense.component_range(0, 0.25).unwrap();
    assert!(prefix[0] <= 0.0 && (prefix[1] - 0.75).abs() < 1e-13);
    assert_eq!(dense.component_range(0, 0.0).unwrap(), [0.0; 2]);
    assert_eq!(dense.component_range(1, 1.0).unwrap(), [-3.0; 2]);
    assert!(dense.component_range(2, 0.5).is_err());
    assert!(dense.component_range(0, -0.01).is_err());
    assert!(dense.component_range(0, 1.01).is_err());
    assert!(dense.component_range(0, f64::NAN).is_err());
}

#[test]
fn observer_receives_only_event_prefix_without_changing_the_original_trajectory() {
    let calls = Cell::new(0usize);
    let rhs = |t: f64, _: &[f64], out: &mut [f64]| {
        calls.set(calls.get() + 1);
        out[0] = 4.0 - 8.0 * t;
    };
    let guard = |_: u64, t: f64, _: &[f64]| t - 0.25;
    let specs = [EventSpec {
        id: 7,
        direction: EventDirection::Rising,
        initial: InitialEvent::Ignore,
    }];
    let mut original = AdaptiveState::new(0.0, &[0.0], 1.0);
    let mut observed = original.clone();
    let baseline = rk45_until_any_event(
        &mut original,
        &rhs,
        &guard,
        &specs,
        1.0,
        1e-11,
        1e-12,
        &PiController::default(),
        4,
        &options(),
        &mut || false,
    )
    .unwrap();
    let original_calls = calls.replace(0);
    let mut seen = Vec::new();
    let result = rk45_until_any_event_observed(
        &mut observed,
        &rhs,
        &guard,
        &specs,
        1.0,
        1e-11,
        1e-12,
        &PiController::default(),
        4,
        &options(),
        &mut || false,
        &mut |dense, accepted_end| {
            seen.push((accepted_end, dense.component_range(0, accepted_end)?));
            Ok(())
        },
    )
    .unwrap();
    assert!(matches!(baseline, EventSetAdvance::Event(_)));
    assert!(matches!(result, EventSetAdvance::Event(_)));
    assert_eq!(seen.len(), 1);
    assert_eq!(seen[0].0, 0.25);
    assert!(
        seen[0].1[1] < 0.751,
        "discarded later peak must not enter the range"
    );
    assert_eq!(observed.t.to_bits(), original.t.to_bits());
    assert_eq!(observed.h.to_bits(), original.h.to_bits());
    assert_eq!(observed.err_prev.to_bits(), original.err_prev.to_bits());
    assert_eq!(observed.u, original.u);
    assert_eq!(
        (observed.accepted, observed.rejected),
        (original.accepted, original.rejected)
    );
    assert_eq!(calls.get(), original_calls);
}

#[test]
fn late_observer_refusal_preserves_the_last_accepted_state() {
    let mut state = AdaptiveState::new(0.0, &[0.0], 0.5);
    let mut opts = options();
    opts.max_step = 0.5;
    let mut visits = 0;
    let result = rk45_until_any_event_observed(
        &mut state,
        &|_, _, out| out[0] = 1.0,
        &|_, _, _| 1.0,
        &[],
        1.0,
        1e-10,
        1e-12,
        &PiController::default(),
        4,
        &opts,
        &mut || false,
        &mut |_, end| {
            visits += 1;
            if end > 0.5 {
                Err(EventError::InvalidSample)
            } else {
                Ok(())
            }
        },
    );
    assert!(result.is_err());
    assert_eq!(visits, 2);
    assert_eq!(state.t, 0.5);
    assert!((state.u[0] - 0.5).abs() < 1e-14);
    assert_eq!(state.accepted, 1);
}

#[test]
fn cancellation_after_observation_does_not_commit_that_step() {
    let mut state = AdaptiveState::new(0.0, &[1.0], 0.5);
    let cancel = Cell::new(false);
    let result = rk45_until_any_event_observed(
        &mut state,
        &|_, _, out| out[0] = 1.0,
        &|_, _, _| 1.0,
        &[],
        0.5,
        1e-10,
        1e-12,
        &PiController::default(),
        4,
        &options(),
        &mut || cancel.get(),
        &mut |_, _| {
            cancel.set(true);
            Ok(())
        },
    )
    .unwrap();
    assert!(
        matches!(result, EventSetAdvance::Stopped(report) if report.status == AdaptiveStatus::Cancelled)
    );
    assert_eq!(state.t, 0.0);
    assert_eq!(state.u, [1.0]);
    assert_eq!(state.accepted, 0);
}
