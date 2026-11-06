//! G0/G1/G3/G4/G5: composition, actual failed-attempt work, and prefix honesty.

use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, BudgetRefusal, CancelGate, ExecMode, StreamKey};

fn with_cx<R>(budget: Budget, f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 38,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            budget,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn unit_rhs(
    _: Interval,
    _: &[Interval],
    out: &mut [Interval],
    _: &Cx<'_>,
) -> Result<(), ValidatedError> {
    out.fill(Interval::point(1.0));
    Ok(())
}

fn quarters() -> TrajectoryOptions {
    TrajectoryOptions {
        initial_step: 0.25,
        ..TrajectoryOptions::default()
    }
}

#[test]
fn g1_constant_flow_composes_uncertain_endpoints_verbatim_and_replays() {
    with_cx(Budget::INFINITE, |cx, _| {
        let initial = [Interval::new(-2.0, -1.0), Interval::new(3.0, 4.0)];
        let run = || {
            validated_trajectory(
                cx,
                2.0,
                3.0,
                &initial,
                &mut |_, _, out, _| {
                    out[0] = Interval::point(2.0);
                    out[1] = Interval::point(-3.0);
                    Ok(())
                },
                &quarters(),
            )
            .unwrap()
        };
        let tube = run();
        assert_eq!(tube, run());
        assert_eq!(tube.interval(), [2.0, 3.0]);
        assert_eq!(tube.initial(), initial);
        assert_eq!(tube.steps().len(), 4);
        for pair in tube.steps().windows(2) {
            assert_eq!(
                pair[0].interval()[1].to_bits(),
                pair[1].interval()[0].to_bits()
            );
            assert_eq!(pair[0].endpoint(), pair[1].initial());
        }
        for step in tube.steps() {
            let t = step.interval()[1] - 2.0;
            assert!(step.endpoint()[0].encloses(Interval::new(-2.0 + 2.0 * t, -1.0 + 2.0 * t)));
            assert!(step.endpoint()[1].encloses(Interval::new(3.0 - 3.0 * t, 4.0 - 3.0 * t)));
        }
        assert_eq!(tube.stats().attempts, tube.steps().len());
        assert_eq!(
            tube.stats().rhs_evaluations,
            tube.steps()
                .iter()
                .map(ValidatedStep::rhs_evaluations)
                .sum::<usize>()
        );
        assert_eq!(
            tube.stats().work,
            3 * initial.len() + tube.steps().iter().map(ValidatedStep::work).sum::<usize>()
        );
    });
}

#[test]
fn g1_inclusion_failure_subdivides_and_encloses_exponential_flow() {
    with_cx(Budget::INFINITE, |cx, _| {
        let options = TrajectoryOptions {
            initial_step: 1.0,
            ..TrajectoryOptions::default()
        };
        let mut calls = 0;
        let tube = validated_trajectory(
            cx,
            0.0,
            1.0,
            &[Interval::point(1.0)],
            &mut |_, state, out, _| {
                calls += 1;
                out[0] = state[0];
                Ok(())
            },
            &options,
        )
        .unwrap();
        assert_eq!(tube.interval(), [0.0, 1.0]);
        assert!(tube.stats().attempts > tube.steps().len());
        assert_eq!(tube.stats().rhs_evaluations, calls);
        assert!(
            calls
                > tube
                    .steps()
                    .iter()
                    .map(ValidatedStep::rhs_evaluations)
                    .sum::<usize>()
        );
        for step in tube.steps() {
            assert!(step.endpoint()[0].contains(step.interval()[1].exp()));
            let midpoint = step.interval()[0].midpoint(step.interval()[1]);
            assert!(step.range(Interval::point(midpoint), cx).unwrap()[0].contains(midpoint.exp()));
        }
    });
}

#[test]
fn g0_width_retries_are_charged_and_minimum_step_is_explicit() {
    with_cx(Budget::INFINITE, |cx, _| {
        let options = TrajectoryOptions {
            initial_step: 0.5,
            step: ValidatedOptions {
                max_width: 0.3,
                ..ValidatedOptions::default()
            },
            ..TrajectoryOptions::default()
        };
        let mut calls = 0;
        let tube = validated_trajectory(
            cx,
            0.0,
            1.0,
            &[Interval::point(0.0)],
            &mut |t, state, out, cx| {
                calls += 1;
                unit_rhs(t, state, out, cx)
            },
            &options,
        )
        .unwrap();
        assert_eq!(tube.interval(), [0.0, 1.0]);
        assert!(tube.endpoint()[0].contains(1.0));
        assert!(tube.stats().attempts > tube.steps().len());
        assert_eq!(tube.stats().rhs_evaluations, calls);
        assert!(
            tube.stats().work > 3 + tube.steps().iter().map(ValidatedStep::work).sum::<usize>()
        );
        assert!(
            tube.steps()
                .iter()
                .all(|s| s.interval()[1] - s.interval()[0] <= 0.125)
        );

        let minimum = TrajectoryOptions {
            min_step: 0.5,
            ..options
        };
        let failure = validated_trajectory(
            cx,
            0.0,
            1.0,
            &[Interval::point(0.0)],
            &mut unit_rhs,
            &minimum,
        )
        .unwrap_err();
        assert_eq!(failure.prefix().unwrap().interval(), [0.0, 0.0]);
        assert_eq!(failure.stats().attempts, 1);
        assert!(matches!(
            failure.stop(),
            TrajectoryStop::MinimumStep {
                last_refusal: ValidatedError::WidthLimit { .. },
                ..
            }
        ));
    });
}

#[test]
fn g0_callback_errors_are_never_reinterpreted_as_adaptive_retries() {
    with_cx(Budget::INFINITE, |cx, _| {
        for error in [
            ValidatedError::RhsRefused("material law is outside its domain"),
            ValidatedError::NoInclusion {
                iterations: 7,
                rhs_evaluations: 9,
            },
            ValidatedError::WidthLimit {
                component: 0,
                width: 9.0,
                limit: 1.0,
            },
            ValidatedError::WorkLimit { limit: 1 },
            ValidatedError::RhsEvaluationLimit { limit: 1 },
            ValidatedError::Arithmetic {
                phase: "callback",
                component: 0,
            },
        ] {
            let mut calls = 0;
            let failure = validated_trajectory(
                cx,
                0.0,
                1.0,
                &[Interval::point(0.0)],
                &mut |_, _, _, _| {
                    calls += 1;
                    Err(error.clone())
                },
                &quarters(),
            )
            .unwrap_err();
            assert_eq!(failure.stop(), &TrajectoryStop::Step(error));
            assert_eq!(failure.stats().attempts, 1);
            assert_eq!(failure.stats().rhs_evaluations, 1);
            assert_eq!(calls, 1);
            assert!(failure.prefix().unwrap().steps().is_empty());
        }
    });
}

#[test]
fn g4_late_refusal_keeps_only_the_identical_completed_prefix() {
    with_cx(Budget::INFINITE, |cx, _| {
        let initial = [Interval::point(0.0)];
        let completed =
            validated_trajectory(cx, 0.0, 0.5, &initial, &mut unit_rhs, &quarters()).unwrap();
        let failure = validated_trajectory(
            cx,
            0.0,
            1.0,
            &initial,
            &mut |time, state, out, cx| {
                if time.lo() >= 0.5 {
                    return Err(ValidatedError::RhsRefused("stop here"));
                }
                unit_rhs(time, state, out, cx)
            },
            &quarters(),
        )
        .unwrap_err();
        let prefix = failure.prefix().unwrap();
        assert_eq!(prefix.interval(), completed.interval());
        assert_eq!(prefix.steps(), completed.steps());
        assert_eq!(prefix.endpoint(), completed.endpoint());
        assert_eq!(failure.stats().attempts, 3);
        assert_eq!(
            failure.stats().rhs_evaluations,
            completed.stats().rhs_evaluations + 1
        );
        assert!(failure.stats().work > completed.stats().work);
        assert_eq!(prefix.stats(), failure.stats());
        assert_eq!(initial, [Interval::point(0.0)]);
    });
}

#[test]
fn g4_cancellation_preserves_the_prefix_and_stops_after_the_current_callback() {
    with_cx(Budget::INFINITE, |cx, gate| {
        gate.request();
        let failure = validated_trajectory(
            cx,
            0.0,
            1.0,
            &[Interval::point(0.0)],
            &mut |_, _, _, _| panic!("pre-cancelled integration evaluated its RHS"),
            &quarters(),
        )
        .unwrap_err();
        assert!(failure.prefix().is_none());
        assert_eq!(
            failure.stop(),
            &TrajectoryStop::Step(ValidatedError::Cancelled)
        );
        assert_eq!(failure.stats(), TrajectoryStats::default());
    });
    with_cx(Budget::INFINITE, |cx, gate| {
        let mut calls = 0;
        let failure = validated_trajectory(
            cx,
            0.0,
            1.0,
            &[Interval::point(0.0)],
            &mut |time, state, out, cx| {
                calls += 1;
                if time.lo() >= 0.5 {
                    gate.request();
                }
                unit_rhs(time, state, out, cx)
            },
            &quarters(),
        )
        .unwrap_err();
        assert_eq!(
            failure.stop(),
            &TrajectoryStop::Step(ValidatedError::Cancelled)
        );
        assert_eq!(failure.prefix().unwrap().interval(), [0.0, 0.5]);
        assert_eq!(failure.stats().attempts, 3);
        assert_eq!(failure.stats().rhs_evaluations, calls);
        assert_eq!(calls, 5);
    });
}

#[test]
fn g0_whole_trajectory_limits_include_refused_attempts() {
    with_cx(Budget::INFINITE, |cx, _| {
        let initial = [Interval::point(0.0)];
        let limited = TrajectoryOptions {
            max_steps: 2,
            ..quarters()
        };
        let failure =
            validated_trajectory(cx, 0.0, 1.0, &initial, &mut unit_rhs, &limited).unwrap_err();
        assert_eq!(failure.stop(), &TrajectoryStop::StepLimit { limit: 2 });
        assert_eq!(failure.prefix().unwrap().interval(), [0.0, 0.5]);
        assert_eq!(failure.stats().attempts, 2);

        let limited = TrajectoryOptions {
            max_attempts: 1,
            step: ValidatedOptions {
                max_width: 0.01,
                ..ValidatedOptions::default()
            },
            ..quarters()
        };
        let failure =
            validated_trajectory(cx, 0.0, 1.0, &initial, &mut unit_rhs, &limited).unwrap_err();
        assert_eq!(failure.stop(), &TrajectoryStop::AttemptLimit { limit: 1 });
        assert_eq!(failure.stats().attempts, 1);
        assert!(failure.prefix().unwrap().steps().is_empty());

        let limited = TrajectoryOptions {
            max_rhs_evaluations: 3,
            ..quarters()
        };
        let mut calls = 0;
        let failure = validated_trajectory(
            cx,
            0.0,
            1.0,
            &initial,
            &mut |t, x, out, cx| {
                calls += 1;
                unit_rhs(t, x, out, cx)
            },
            &limited,
        )
        .unwrap_err();
        assert_eq!(
            failure.stop(),
            &TrajectoryStop::RhsEvaluationLimit { limit: 3 }
        );
        assert_eq!(failure.prefix().unwrap().interval(), [0.0, 0.25]);
        assert_eq!(failure.stats().rhs_evaluations, 3);
        assert_eq!(calls, 3);

        let one =
            validated_trajectory(cx, 0.0, 0.25, &initial, &mut unit_rhs, &quarters()).unwrap();
        let limited = TrajectoryOptions {
            max_work: one.stats().work + 6,
            ..quarters()
        };
        let failure =
            validated_trajectory(cx, 0.0, 1.0, &initial, &mut unit_rhs, &limited).unwrap_err();
        assert_eq!(
            failure.stop(),
            &TrajectoryStop::WorkLimit {
                limit: limited.max_work
            }
        );
        assert_eq!(failure.prefix().unwrap().interval(), [0.0, 0.25]);
        assert!(failure.stats().work > one.stats().work);
        assert!(failure.stats().work <= limited.max_work);
    });
}

#[test]
fn g4_ambient_budget_is_shared_across_steps_and_not_readmitted_per_attempt() {
    with_cx(
        Budget {
            poll_quota: 40,
            ..Budget::INFINITE
        },
        |cx, _| {
            let failure = validated_trajectory(
                cx,
                0.0,
                1.0,
                &[Interval::point(0.0)],
                &mut unit_rhs,
                &quarters(),
            )
            .unwrap_err();
            assert!(matches!(
                failure.stop(),
                TrajectoryStop::Step(ValidatedError::Budget(BudgetRefusal::PollsExhausted {
                    quota: 40,
                    ..
                }))
            ));
            assert_eq!(failure.prefix().unwrap().interval(), [0.0, 0.25]);
        },
    );
    with_cx(
        Budget {
            cost_quota: Some(1),
            ..Budget::INFINITE
        },
        |cx, _| {
            let failure = validated_trajectory(
                cx,
                0.0,
                1.0,
                &[Interval::point(0.0)],
                &mut |_, _, _, _| panic!("failed admission evaluated the RHS"),
                &quarters(),
            )
            .unwrap_err();
            assert!(failure.prefix().is_none());
            assert!(matches!(
                failure.stop(),
                TrajectoryStop::Step(ValidatedError::Budget(
                    BudgetRefusal::CostPlanExceedsQuota { .. }
                ))
            ));
        },
    );
}

#[test]
fn g0_zero_horizon_and_subminimum_final_remainder_have_exact_coverage() {
    with_cx(Budget::INFINITE, |cx, _| {
        let initial = [Interval::new(-0.5, 0.5)];
        let tube = validated_trajectory(
            cx,
            3.0,
            3.0,
            &initial,
            &mut |_, _, _, _| panic!("zero horizon evaluated the RHS"),
            &quarters(),
        )
        .unwrap();
        assert_eq!(tube.interval(), [3.0, 3.0]);
        assert_eq!(tube.initial(), initial);
        assert_eq!(tube.endpoint(), initial);
        assert!(tube.steps().is_empty());
        assert_eq!(
            tube.stats(),
            TrajectoryStats {
                attempts: 0,
                rhs_evaluations: 0,
                work: 3
            }
        );

        let options = TrajectoryOptions {
            min_step: 0.25,
            ..quarters()
        };
        let tube = validated_trajectory(
            cx,
            0.0,
            0.3,
            &[Interval::point(0.0)],
            &mut unit_rhs,
            &options,
        )
        .unwrap();
        assert_eq!(tube.interval(), [0.0, 0.3]);
        assert_eq!(tube.steps().len(), 2);
        assert_eq!(tube.steps()[1].interval(), [0.25, 0.3]);
        assert!(tube.endpoint()[0].contains(0.3));
    });
}

#[test]
fn g0_invalid_inputs_time_resolution_and_arithmetic_never_publish_a_bad_step() {
    with_cx(Budget::INFINITE, |cx, _| {
        for options in [
            TrajectoryOptions {
                min_step: 0.0,
                ..quarters()
            },
            TrajectoryOptions {
                initial_step: f64::INFINITY,
                ..quarters()
            },
            TrajectoryOptions {
                max_attempts: 0,
                ..quarters()
            },
            TrajectoryOptions {
                max_steps: usize::MAX,
                ..quarters()
            },
        ] {
            let failure = validated_trajectory(
                cx,
                0.0,
                1.0,
                &[Interval::point(0.0)],
                &mut |_, _, _, _| panic!("invalid options called RHS"),
                &options,
            )
            .unwrap_err();
            assert!(failure.prefix().is_none());
            assert!(matches!(
                failure.stop(),
                TrajectoryStop::Step(ValidatedError::InvalidOptions(_))
            ));
        }
        let failure = validated_trajectory(
            cx,
            1.0e16,
            1.0e16 + 4.0,
            &[Interval::point(0.0)],
            &mut unit_rhs,
            &quarters(),
        )
        .unwrap_err();
        assert!(matches!(
            failure.stop(),
            TrajectoryStop::TimeResolution { .. }
        ));
        assert!(failure.prefix().unwrap().steps().is_empty());
        assert_eq!(failure.stats().attempts, 0);

        let failure = validated_trajectory(
            cx,
            0.0,
            1.0,
            &[Interval::point(f64::MAX)],
            &mut |_, _, _, _| panic!("overflowing inflation called RHS"),
            &quarters(),
        )
        .unwrap_err();
        assert!(matches!(
            failure.stop(),
            TrajectoryStop::Step(ValidatedError::Arithmetic { .. })
        ));
        assert_eq!(failure.stats().attempts, 1);
        assert!(failure.prefix().unwrap().steps().is_empty());
        assert_eq!(failure.stats().rhs_evaluations, 0);
    });
}
