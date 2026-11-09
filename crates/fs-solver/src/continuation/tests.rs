use super::*;
use std::cell::Cell;

fn controls() -> ContinuationControls {
    ContinuationControls {
        initial_increment: 0.5,
        minimum_increment: 1.0 / 128.0,
        maximum_increment: 0.5,
        maximum_attempts: 64,
    }
}

#[test]
fn nonlinear_branch_keeps_only_admitted_seeds_and_reaches_exact_target() {
    // The corrector's independently known branch is q^3 + q = parameter.
    // Its local convergence basin is deliberately bounded; the scheduler must
    // recover without altering either that root equation or an accepted seed.
    let mut targets = Vec::new();
    let answer = follow(
        (0.0, 0.0_f64),
        controls(),
        |target, &(before, q)| {
            assert!((q * q * q + q - before).abs() < 1e-14);
            targets.push(target);
            if target - before > 0.125 {
                return Err(CorrectorRefusal::Retry("outside local basin"));
            }
            let (mut lo, mut hi) = (q, 1.0);
            for _ in 0..60 {
                let mid = 0.5 * (lo + hi);
                if mid * mid * mid + mid < target {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let next = 0.5 * (lo + hi);
            assert!((next * next * next + next - target).abs() < 1e-14);
            Ok((target, next))
        },
        || Ok::<_, &str>(()),
    )
    .unwrap();
    assert_eq!(answer.state.0, 1.0);
    assert!((answer.state.1 - 0.6823278038280193).abs() < 1e-14);
    assert_eq!(answer.report.accepted_steps, 8);
    assert!(answer.report.rejected_steps > 0);
    assert_eq!(answer.report.attempts, targets.len());
    assert_eq!(targets.last(), Some(&1.0));
    assert_eq!(answer.report.smallest_accepted_increment, 0.125);
}

#[test]
fn original_source_error_stops_immediately_and_never_becomes_a_retry() {
    let calls = Cell::new(0);
    let error = follow(
        (),
        controls(),
        |_, _| {
            calls.set(calls.get() + 1);
            Err::<(), _>(CorrectorRefusal::Stop("material domain"))
        },
        || Ok::<_, &str>(()),
    )
    .unwrap_err();
    assert_eq!(calls.get(), 1);
    assert_eq!(
        error,
        ContinuationError::Source {
            parameter: 0.5,
            source: "material domain"
        }
    );
}

#[test]
fn unresolved_fold_preserves_its_original_refusal_and_never_returns_a_partial_root() {
    let error = follow(
        0.0,
        controls(),
        |target, _| {
            if target <= 0.5 {
                Ok(target)
            } else {
                Err(CorrectorRefusal::Retry("no equilibrium beyond fold"))
            }
        },
        || Ok::<_, &str>(()),
    )
    .unwrap_err();
    assert!(matches!(
        error,
        ContinuationError::IncrementLimit {
            accepted_parameter: 0.5,
            source: "no equilibrium beyond fold",
            ..
        }
    ));
}

#[test]
fn final_checkpoint_cancels_before_publication_and_retry_is_repeatable() {
    let ready = Cell::new(false);
    let one = ContinuationControls {
        initial_increment: 1.0,
        maximum_increment: 1.0,
        maximum_attempts: 1,
        ..controls()
    };
    let error = follow(
        7,
        one,
        |p, &old| {
            assert_eq!(p, 1.0);
            assert_eq!(old, 7);
            ready.set(true);
            Ok(9)
        },
        || {
            if ready.get() {
                Err("cancelled")
            } else {
                Ok(())
            }
        },
    )
    .unwrap_err();
    assert_eq!(
        error,
        ContinuationError::Source {
            parameter: 1.0,
            source: "cancelled"
        }
    );
    let again = || {
        follow(
            7,
            one,
            |_, &old| Ok::<_, CorrectorRefusal<&str>>(old + 2),
            || Ok(()),
        )
        .unwrap()
    };
    assert_eq!(again(), again());
    assert_eq!(again().state, 9);
}

#[test]
fn attempt_cap_and_small_final_remainder_have_exact_boundary_semantics() {
    let mut c = controls();
    c.initial_increment = 0.3;
    c.maximum_increment = 0.3;
    c.minimum_increment = 0.2;
    c.maximum_attempts = 4;
    let solve = |p, _: &f64| Ok::<_, CorrectorRefusal<&str>>(p);
    let result = follow(0.0, c, solve, || Ok(())).unwrap();
    assert_eq!(result.state, 1.0);
    assert_eq!(result.report.attempts, 4);
    assert!(result.report.smallest_accepted_increment < c.minimum_increment);
    c.maximum_attempts = 3;
    assert!(matches!(
        follow(0.0, c, solve, || Ok(())),
        Err(ContinuationError::AttemptLimit { attempts: 3, .. })
    ));
}

#[test]
fn invalid_controls_make_no_corrector_call() {
    let good = controls();
    for c in [
        ContinuationControls {
            minimum_increment: 0.0,
            ..good
        },
        ContinuationControls {
            initial_increment: f64::NAN,
            ..good
        },
        ContinuationControls {
            maximum_increment: 1.1,
            ..good
        },
        ContinuationControls {
            maximum_attempts: 0,
            ..good
        },
        ContinuationControls {
            maximum_attempts: 257,
            ..good
        },
    ] {
        let error = follow(
            (),
            c,
            |_, _| -> Result<(), CorrectorRefusal<&str>> {
                panic!("invalid controls reached the source")
            },
            || Ok(()),
        )
        .unwrap_err();
        assert!(matches!(error, ContinuationError::Invalid(_)));
    }
}
